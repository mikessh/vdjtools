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

class ImSeqParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected ImSeqParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.ImSeq, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        def splitString = clonotypeString.split("[: \t]")
        def count = splitString[3].toInteger()
        def freq = 0

        def cdr3nt = splitString[1].toUpperCase()
        def cdr3aa = toUnifiedCdr3Aa(translate(cdr3nt))

        String v, d, j
        (v, j, d) = extractVDJ([splitString[0], splitString[2], "."])


        def segmPoints = [-1, -1, -1, -1] as int[]

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa), isComplete = true

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)
    }
}
