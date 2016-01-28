/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */


package com.antigenomics.vdjtools.io.parser

import com.antigenomics.vdjtools.misc.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample

import static com.antigenomics.vdjtools.misc.CommonUtil.*

/**
 * A clonotype parser implementation that handles Adaptive Biotechnologies (tm) immunoSEQ (tm) assay
 * output format, see
 * {@url http://www.adaptivebiotech.com/content/immunoseq-0}
 */
class ImmunoSeqParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected ImmunoSeqParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.ImmunoSeq, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        /*
             $id   | content                    | description
            -------|----------------------------|------------------------------------------------------------------------
             00    | sequencing read            | this is raw sequence of a read. It could contain esither full or partial (5') CDR3 sequence
             01    | CDR3 amino acid sequence   | standard AA sequence or blank, in case read was too short to cover CDR3
             02    | Count                      | default
             03    | Fraction                   | fraction that also accounts for incomplete reads
             04-05 | unused                     |
             06    | V segments                 | V "family", non-conventional naming
             07-12 | unused                     |
             13    | D segments                 | D "family", non-conventional naming
             14-19 | unused                     |
             20    | J segments                 | J "family", non-conventional naming
             21-31 | unused                     |
             32    | CDR3 start                 | used to extract CDR3 nucleotide sequence from $00
             33    | V end                      |
             34    | D start                    |
             35    | D end                      |
             36    | J start                    |
             37    | unused                     |
             38+   | unused                     |
             
             - note that there is no "J end" (perhaps due to J segment identification issues), so
              here we use $32 + length($01) * 3
          */

        def splitString = clonotypeString.split("\t")

        def count = splitString[2].toInteger()
        def freq = splitString[3].toDouble()

        // This field is used to extract CDR3 region, as the data contains unprocessed sequences in $00 field.
        // For data that was already processed by VDJtools, an extracted CDR3 sequence is stored to $00 field
        // and $32 is set to "." to indicate that no additional CDR3 nucleotide sequence extraction is needed.
        def cdr3start = splitString[32].isInteger() ?
                splitString[32].toInteger() :
                -1

        def cdr3nt = splitString[0]
        def cdr3aa = splitString[1]

        def jStart = splitString[36].toInteger()

        boolean isComplete

        if (cdr3start >= 0) {
            if (cdr3aa.length() > 0) {
                int to = cdr3start + cdr3aa.length() * 3
                // see https://github.com/mikessh/vdjtools/issues/30 for the reason for workaround
                isComplete = to <= cdr3nt.length()
                cdr3nt = isComplete ? cdr3nt.substring(cdr3start, to) : cdr3nt.substring(cdr3start)// in-frame
            } else {
                // it seems to be hard to get conventional out-of-frame translation here
                // but we'll try to reconstruct it
                if (jStart >= 0) {
                    def jRef = getJReferencePoint(cdr3nt.substring(jStart))
                    if (jRef >= 0) {
                        cdr3nt = cdr3nt.substring(cdr3start, jStart + jRef + 4)
                        cdr3aa = translate(cdr3nt).replaceAll(/([atgc#\?])+/, "~") // to unified look
                    }
                }
                isComplete = cdr3aa.length() > 0
            }
        } else {
            isComplete = false
        }

        String v, d, j
        (v, d, j) = extractVDJ(splitString[[7, 14, 21]]).collect {
            it.toLowerCase() == "unresolved" ? "." : it
        }

        boolean inFrame = cdr3aa.length() > 0 && inFrame(cdr3aa),
                noStop = noStop(cdr3aa)

        // Correctly record segment points
        cdr3start = cdr3start < 0 ? 0 : cdr3start

        def segmPoints = [
                splitString[33].toInteger() - 1 - cdr3start,
                splitString[34].toInteger() - cdr3start,
                splitString[35].toInteger() - 1 - cdr3start,
                jStart - cdr3start].collect { it < 0 ? -1 : it } as int[]

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)
    }
}
