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
import com.antigenomics.vdjtools.util.CommonUtil

class Clonotype implements Countable {
    VSegmentTable parentVSegmentTable = null
    final Set<Mutation> mutations = new HashSet<Mutation>(),
                        alleles = new HashSet<Mutation>(),
                        shms = new HashSet<Mutation>()

    int count
    double freq

    final String v, d, j
    final String cdr1nt, cdr2nt, cdr3nt,
                 cdr1aa, cdr2aa, cdr3aa

    final String blank = "."

    private String key

    String getKey() {
        key
    }

    final boolean inFrame, isComplete, noStop

    Clonotype changeCount(int newCount, long total) {
        new Clonotype(newCount, newCount / (double) total,
                v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, isComplete, noStop)
    }

    private Clonotype(int count, double freq,
                      String v, String d, String j,
                      String cdr1nt, String cdr2nt, String cdr3nt,
                      String cdr1aa, String cdr2aa, String cdr3aa,
                      boolean inFrame, boolean isComplete, boolean noStop) {
        this.count = count
        this.freq = freq
        this.v = v
        this.d = d
        this.j = j
        this.cdr1nt = cdr1nt
        this.cdr2nt = cdr2nt
        this.cdr3nt = cdr3nt
        this.cdr1aa = cdr1aa
        this.cdr2aa = cdr2aa
        this.cdr3aa = cdr3aa
        this.key = key
        this.inFrame = inFrame
        this.isComplete = isComplete
        this.noStop = noStop
    }

    static Clonotype parseClonotype(String clonotypeString, Software software) {
        try {
            switch (software) {
                case Software.MiTcr:
                    return parseMiTcrClonotype(clonotypeString)
                case Software.IgBlast:
                    return parseIgBlastClonotype(clonotypeString)
                case Software.Simple:
                    return parseSimpleClonotype(clonotypeString)
            }
        } catch (Exception e) {
            throw new Exception("Error parsing clonotype from $clonotypeString")
        }
        throw new UnsupportedOperationException("Don't know how to parse $software data")
    }

    static Clonotype parseSimpleClonotype(String clonotypeString) {
        def splitString = clonotypeString.split("\t")

        def count = splitString[0].toInteger()
        def freq = splitString[1].toDouble()

        String cdr1nt = null, cdr2nt = null, cdr3nt, cdr1aa = null, cdr2aa = null, cdr3aa
        cdr3nt = splitString[2]
        cdr3aa = splitString[3]

        String v, d, j
        (v, d, j) = extractVDJ(splitString[4..-1])

        boolean inFrame = !(cdr3aa.contains("~") || cdr3aa.contains("?")),
                noStop = !cdr3aa.contains("*"), isComplete = true

        def clonotype = new Clonotype(count, freq,
                v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete)

        clonotype.key = [v, cdr3nt].join("_")

        clonotype
    }

    static Clonotype parseMiTcrClonotype(String clonotypeString) {
        def splitString = clonotypeString.split("\t")

        def count = splitString[0].toInteger()
        def freq = splitString[1].toDouble()

        String cdr1nt = null, cdr2nt = null, cdr3nt, cdr1aa = null, cdr2aa = null, cdr3aa
        cdr3nt = splitString[2]
        cdr3aa = splitString[5]

        String v, d, j
        (v, d, j) = extractVDJ(splitString[[7, 11, 9]])

        boolean inFrame = !cdr3aa.contains("~"), noStop = !cdr3aa.contains("*"), isComplete = true

        def clonotype = new Clonotype(count, freq,
                v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete)

        clonotype.key = [v, cdr3nt].join("_")

        clonotype
    }

    static Clonotype parseMiGecClonotype(String clonotypeString) {
        /*
            0  Count
            1  Percentage
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

        String v, d, j
        (v, d, j) = extractVDJ(splitString[4..6])

        boolean inFrame = !cdr3aa.contains("~"), noStop = !cdr3aa.contains("*"), isComplete = true

        def clonotype = new Clonotype(count, freq,
                v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete)

        clonotype.key = [v, cdr3nt].join("_")

        clonotype
    }

    static Clonotype parseIgBlastClonotype(String clonotypeString) {
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
        (v, d, j) = extractVDJ(splitString[13..15])

        boolean inFrame, noStop, isComplete
        (inFrame, noStop, isComplete) = splitString[10..12].collect { it.toBoolean() }

        def clonotype = new Clonotype(count, freq,
                v, d, j,
                cdr1nt, cdr2nt, cdr3nt,
                cdr1aa, cdr2aa, cdr3aa,
                inFrame, noStop, isComplete)

        if (splitString[19] != ".")
            clonotype.mutations.addAll(splitString[19].split("\\|").collect { String mutString ->
                Mutation.parseIgBlastMutation(mutString, clonotype)
            })

        def key = [v, cdr3nt,
                   clonotype.mutations.collect {
                       Mutation mpd -> "$mpd.ntPos:$mpd.fromNt>$mpd.toNt"
                   }.join("|")].join("_")

        clonotype.key = key

        clonotype
    }

    /**
     * Internal util, extracts most probable allele
     */
    private static List<String> extractVDJ(List<String> vdj) {
        vdj.collect {
            def major = it.split(",")[0]
            major = major.split("\\*")[0] // trim allele if present
            major = major.replaceAll("\"", "")  // zap characters introduced by opening file in Excel
            major.length() > 0 ? major : ""
        }
    }

    String getSubSequence(int from, int to) {
        CommonUtil.getSubSequence(cdr3nt, from, to)
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

    final static String HEADER =
            "display_name\tcount\tfreq\t" +
                    "cdr1nt\tcdr2nt\tcdr3nt\t" +
                    "cdr1aa\tcdr2aa\tcdr3aa\t" +
                    "inFrame\tnoStop\tisComplete\t" +
                    "V\tD\tJ\t" +
                    "shms\talleles"
    final static String NODE_HEADER = "key\t" + HEADER

    String nodeString() {
        key + "\t" + toString()
    }

    @Override
    String toString() {
        [displayName, count, freq,
         cdr1nt, cdr2nt, cdr3nt,
         cdr1aa, cdr2aa, cdr3aa,
         inFrame, noStop, isComplete,
         v, d, j,
         shms.collect { it.key }.join("|"), alleles.collect { it.key }.join("|")].join("\t")
    }

    //
    // Table output
    //

    void print(PrintWriter printWriter, Software software) {
        printWriter.println(software.printFields.collect { this."$it" }.join("\t"))
    }
}
