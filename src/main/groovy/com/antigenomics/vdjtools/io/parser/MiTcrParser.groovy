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

/**
 * A clonotype parser implementation that handles output from MiTCR software, see
 * {@url http://mitcr.milaboratory.com/}
 */
public class MiTcrParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    protected MiTcrParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.MiTcr, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        def splitString = clonotypeString.split(software.delimiter)

        def count = splitString[0].toInteger()
        def freq = splitString[1].toDouble()

        def cdr3nt = splitString[2]
        def cdr3aa = toUnifiedCdr3Aa(splitString[5]) // replace ~

        String v, d, j
        (v, d, j) = extractVDJ(splitString[[7, 11, 9]])

        def segmPoints = [splitString[12].toInteger(),
                          splitString[13].isInteger() ? splitString[13].toInteger() : -1,
                          splitString[14].isInteger() ? splitString[14].toInteger() : -1,
                          splitString[15].toInteger()] as int[]

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa),
                isComplete = true

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j, cdr3nt, cdr3aa,
                inFrame, noStop, isComplete)
    }
}
