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
import com.antigenomics.vdjtools.util.CommonUtil

import static com.antigenomics.vdjtools.util.CommonUtil.inFrame
import static com.antigenomics.vdjtools.util.CommonUtil.noStop

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

        def cdr3start = splitString[58].isInteger() ?
                splitString[58].toInteger() :
                -1

        String cdr1nt = null, cdr2nt = null, cdr3nt, cdr1aa = null, cdr2aa = null, cdr3aa

        cdr3nt = splitString[14]
        cdr3aa = CommonUtil.translate(cdr3nt)

        String v, d, j
        (v, d, j) = CommonUtil.extractVDJ(splitString[3..5])


        def segmPoints = [
                splitString[58],
                splitString[76],
                splitString[77],
                splitString[106]].collect { it.isInteger() ? (it.toInteger() - cdr3start) : -1 } as int[]

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa), isComplete = cdr3aa.length() > 0


        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete,
                new HashSet<>())
    }
}
