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
class ImmunoSeqV3Parser extends ClonotypeStreamParser {
    private boolean initialized = false
    private int countColumn, freqColumn, cdr3StartColumn,
                cdr3ntColumn, cdr3aaColumn, cdr3LenColumn,
                jStartColumn, inFrameColumn,
                vColumn, dColumn, jColumn,
                vEndColumn, dStartColumn, dEndColumn

    /**
     * {@inheritDoc}
     */
    protected ImmunoSeqV3Parser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.ImmunoSeq, sample)
    }

    /**
     * Performs parser initialization based on the header string.
     *
     * @throws RuntimeException if header line doesn't contain required columns
     */
    private synchronized void ensureInitialized() {
        if (initialized)
            return

        // Parsing header line to determine positions of certain columns with clones properties

        String headerLine = this.header[0];
        String[] splitHeaderLine = headerLine.split(software.delimiter)

        countColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("reads") }
        freqColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("frequency") }
        cdr3StartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v_index") }
        cdr3LenColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("cdr3_length") }
        cdr3ntColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("rearrangement") }
        inFrameColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("frame_type") }
        cdr3aaColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("amino_acid") }
        vColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("v_family") }
        dColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d_family") }
        jColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j_family") }
        vEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("n1_index") }
        dStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("d_index") }
        dEndColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("n2_index") }
        jStartColumn = splitHeaderLine.findIndexOf { it.equalsIgnoreCase("j_index") }

        if ([countColumn, freqColumn,
             cdr3StartColumn,  cdr3LenColumn,
             cdr3ntColumn, cdr3aaColumn,
             vColumn, dColumn, jColumn,
             vEndColumn, dStartColumn, dEndColumn, jStartColumn,
             inFrameColumn].any { it < 0 })
            throw new RuntimeException("Some mandatory columns are absent in the input file.")

        // Initialized
        initialized = true
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        ensureInitialized()

        def splitString = clonotypeString.split("\t")

        def count = splitString[countColumn].toInteger()
        def freq = splitString[freqColumn].toDouble()

        def cdr3start = splitString[cdr3StartColumn].toInteger(),
            cdr3Len = splitString[cdr3LenColumn].toInteger()

        def inFrame = splitString[inFrameColumn].equalsIgnoreCase("in")

        def cdr3nt = splitString[cdr3ntColumn][cdr3start..<(cdr3start + cdr3Len)]
        def cdr3aa = toUnifiedCdr3Aa(inFrame ? splitString[cdr3aaColumn] : translate(cdr3nt))

        String v, d, j
        (v, d, j) = extractVDJ(splitString[[vColumn, dColumn, jColumn]])

        boolean noStop = noStop(cdr3aa)

        // Correctly record segment points
        cdr3start = cdr3start < 0 ? 0 : cdr3start

        def segmPoints = [
                splitString[vEndColumn].toInteger() - 1 - cdr3start,
                splitString[dStartColumn].toInteger() - cdr3start,
                splitString[dEndColumn].toInteger() - 1 - cdr3start,
                splitString[jStartColumn].toInteger() - cdr3start].collect { it < 0 ? -1 : it } as int[]

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, true)
    }
}
