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
 * A clonotype parser implementation that handles simple tab-delimited input, see
 * {@url https://github.com/mikessh/vdjtools/wiki/Input#simple}
 */
public class SimpleParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    public SimpleParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.Simple, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        def splitString = clonotypeString.split(software.delimiter)

        def count = splitString[0].toInteger()
        def freq = splitString[1].toDouble()

        String cdr1nt = null, cdr2nt = null, cdr3nt, cdr1aa = null, cdr2aa = null, cdr3aa

        cdr3nt = splitString[2]
        cdr3aa = toUnifiedCdr3Aa(splitString[3])

        String v, d, j
        (v, d, j) = extractVDJ(splitString[4..6])

        boolean inFrame = inFrame(cdr3aa),
                noStop = noStop(cdr3aa),
                isComplete = true


        def segmPoints = (7..10).collect {
            (splitString.size() <= it || !splitString[it].isInteger()) ? -1 : splitString[it].toInteger()
        } as int[]

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete,
                null)
    }
}
