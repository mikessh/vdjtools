/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
