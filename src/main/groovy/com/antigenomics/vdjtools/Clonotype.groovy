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

package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.segment.SegmentData
import com.antigenomics.vdjtools.segment.VSegmentTable

class Clonotype {
    VSegmentTable parentVSegmentTable = null
    final Set<Mutation> mutations = new HashSet<Mutation>(),
                        alleles = new HashSet<Mutation>(),
                        shms = new HashSet<Mutation>()

    int count
    double freq

    final String v, d, j
    final String cdr1nt, cdr2nt, cdr3nt, cdr1aa, cdr2aa, cdr3aa, key

    final boolean inFrame, isComplete, noStop

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

    Clonotype(String clonotypeString) {
        def splitString = clonotypeString.split("\t")

        count = splitString[2].toInteger()
        freq = splitString[3].toDouble()

        (cdr1nt, cdr2nt, cdr3nt, cdr1aa, cdr2aa, cdr3aa) = splitString[4..9]

        (v, d, j) = splitString[13..15]

        (inFrame, noStop, isComplete) = splitString[10..12].collect { it.toBoolean() }

        if (splitString[19] != ".")
            mutations.addAll(splitString[19].split("\\|").collect { String mutString ->
                Mutation.parseIgBlastMutation(mutString, this)
            })

        key = [cdr3nt, mutations.collect { Mutation mpd -> "$mpd.ntPos:$mpd.fromNt>$mpd.toNt" }.join("|")].join("_")
    }

    String getSubSequence(int from, int to) {
        Util.getSubSequence(cdr3nt, from, to)
    }

    SegmentData getVSegmentData() {
        parentVSegmentTable.getSegmentData(v)
    }

    //
    // Display
    //

    String getDisplayName() {
        "V" + v[2] + v[4..-4] + ":" + cdr3aa + ":S" + shms.size()
    }

    final static NODE_HEADER =
            "display_name\tcount\tfreq\t" +
                    "cdr1nt\tcdr2nt\tcdr3nt\t" +
                    "cdr1aa\tcdr2aa\tcdr3aa\t" +
                    "inFrame\tnoStop\tisComplete\t" +
                    "V\tD\tJ\t" +
                    "shms\talleles"

    @Override
    String toString() {
        [displayName, count, freq,
         cdr1nt, cdr2nt, cdr3nt,
         cdr1aa, cdr2aa, cdr3aa,
         inFrame, noStop, isComplete,
         v, d, j,
         shms.collect { it.key }.join("|"), alleles.collect { it.key }.join("|")].join("\t")
    }
}
