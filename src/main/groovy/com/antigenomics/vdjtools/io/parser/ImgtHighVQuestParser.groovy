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

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample

import static com.antigenomics.vdjtools.util.CommonUtil.*

class ImgtHighVQuestParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected ImgtHighVQuestParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.ImgtHighVQuest, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        def splitString = clonotypeString.split("\t")

        if (splitString.length < 107)
            return null

        def count = 1
        def freq = 0

        def cdr3start = splitString[60].isInteger() ?
                splitString[60].toInteger() :
                -1 // this is called "junction start" here. Junction = CDR3 + conserved C, F/W

        def cdr3nt = splitString[15].toUpperCase()

        if (!(cdr3nt =~ /^[ATGCatgc]+$/))
            return null // no N's allowed

        def cdr3aa = toUnifiedCdr3Aa(translate(cdr3nt))

        String v, d, j
        (v, j, d) = extractVDJ(splitString[3..5]).collect {
            def splitRecord = it.split(" ")
            splitRecord.length > 1 ? splitRecord[1] : splitRecord[0]
        }

        def segmPoints = [
                splitString[63],
                splitString[76],
                splitString[77],
                splitString[106]
        ].collect {
            it.isInteger() ? (it.toInteger() - cdr3start) : -1 // subtract cdr3start
        }.collect {
            (it >= 0 && it < cdr3nt.length()) ? it : -1 // sometimes segment bounds appear out of junction region
        } as int[]

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa), isComplete = cdr3aa.length() > 0


        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)
    }
}
