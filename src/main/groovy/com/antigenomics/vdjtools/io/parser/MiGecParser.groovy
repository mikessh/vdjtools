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

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.CommonUtil

/**
 * A clonotype parser implementation that handles output from MiGEC software, see
 * {@url https://github.com/mikessh/migec}
 */
public class MiGecParser extends ClonotypeStreamParser {
    /**
     * {@inheritDoc}
     */
    public MiGecParser(Iterator<String> innerIter, Sample sample) {
        super(innerIter, Software.MiGec, sample)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected Clonotype innerParse(String clonotypeString) {
        /*
             0  Count
             1  Fraction
             2  CDR3 nucleotide sequence
             3  CDR3 amino acid sequence
             4  V segments
             5  J segments
             6  D segments
             7  Last V nucleotide position
             8  First D nucleotide position
             9  Last D nucleotide position
             10 First J nucleotide position
             11 Good events
             12 Total events
             13 Good reads
             14 Total reads
          */

        def splitString = clonotypeString.split("\t")

        def count = splitString[0].toInteger()
        def freq = splitString[1].toDouble()

        String cdr1nt = null, cdr2nt = null, cdr3nt, cdr1aa = null, cdr2aa = null, cdr3aa
        cdr3nt = splitString[2]
        cdr3aa = splitString[3]


        String v, j, d
        (v, j, d) = CommonUtil.extractVDJ(splitString[4..6])

        boolean inFrame = !cdr3aa.contains("?"), noStop = !cdr3aa.contains("*"), isComplete = true

        def segmPoints = [
                splitString[7].toInteger(),
                splitString[8].isInteger() ? splitString[8].toInteger() : -1,
                splitString[9].isInteger() ? splitString[9].toInteger() : -1,
                splitString[10].toInteger()] as int[]

        new Clonotype(sample, count, freq,
                segmPoints, v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete,
                new HashSet<>())
    }
}
