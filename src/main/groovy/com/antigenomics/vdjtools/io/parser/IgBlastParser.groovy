package com.antigenomics.vdjtools.io.parser

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Mutation
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.CommonUtil

/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */
class IgBlastParser extends ClonotypeStreamParser {
    IgBlastParser(Iterator<String> innerIter, Software software, Sample sample) {
        super(innerIter, software, sample)
    }

    @Override
    protected Clonotype innerParse(String clonotypeString) {
        /*
        0	1	2	3
        #reads_count	reads_percent	events_count	events_percent
        4	5	6	7	8	9
        cdr1nt	cdr2nt	cdr3nt	cdr1aa	cdr2aa	cdr3aa
        10	11	12
        inFrame	noStop	complete
        13	14	15
        vSegment	dSegment	jSegment
        16	17	18
        cdr1q	cdr2q	cdr3q
        19
        mutations
         */

        def splitString = clonotypeString.split("\t")

        def count = splitString[2].toInteger()
        def freq = splitString[3].toDouble()

        String cdr1nt, cdr2nt, cdr3nt, cdr1aa, cdr2aa, cdr3aa
        (cdr1nt, cdr2nt, cdr3nt, cdr1aa, cdr2aa, cdr3aa) = splitString[4..9]

        String v, d, j
        (v, d, j) = CommonUtil.extractVDJ(splitString[13..15])

        boolean inFrame, noStop, isComplete
        (inFrame, noStop, isComplete) = splitString[10..12].collect { it.toBoolean() }

        def mutations = new HashSet<Mutation>()

        def clonotype = new Clonotype(sample,
                count, freq,
                [-1, -1, -1, -1] as int[], v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete,
                mutations)

        if (splitString[19] != ".")
            mutations.addAll(splitString[19].split("\\|").collect { String mutString ->
                parseIgBlastMutation(mutString, clonotype)
            })

        return clonotype
    }

    private static Mutation parseIgBlastMutation(String mutationString, Clonotype parent) {
        def splitString = mutationString.split(",")

        def ntString = splitString[1]
        def splitNTString = ntString.split("[:>]")
        def ntPos = splitNTString[0].toInteger()
        char fromNt = splitNTString[1]
        char toNt = splitNTString[2]

        def aaString = splitString[2]
        def splitAAString = aaString.split("[:>]")
        def aaPos = splitAAString[0].toInteger()
        char fromAa = splitAAString[1]
        char toAa = splitAAString[2]

        def region = splitString[3]

        new Mutation(region,
                ntPos, aaPos,
                fromAa, toAa, fromNt, toNt,
                true, parent)
    }
}
