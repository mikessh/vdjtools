package com.antigenomics.vdjtools.segment

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Util
import com.antigenomics.vdjtools.io.FastaReader

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
class VSegmentTable {
    final Map<String, FrequencyCounter> countersByName = new HashMap()
    final Map<String, SegmentData> segmentByName = new HashMap<>()

    VSegmentTable(String species) {
        def segmentReader = new FastaReader(Util.resourceStreamReader("segments/${species}_V.fa"))
        def segmentToSequenceMap = segmentReader.collectEntries {
            [(it.header): it.sequence]
        }

        Util.resourceStreamReader("segments/${species}_regions.txt").splitEachLine("\t") { List<String> splitLine ->
            def segmentName = splitLine[0]
            def sequence = segmentToSequenceMap[segmentName]

            def segmentData = new SegmentData(segmentName, sequence,
                    new IntRange(1, 10).step(2).collect { int i ->
                        new Range(splitLine[i].toInteger() - 1, // 1-based
                                splitLine[i + 1].toInteger()    // non-inclusive
                        )
                    }
            )

            countersByName.put(segmentName, new FrequencyCounter())
            segmentByName.put(segmentName, segmentData)
        }
    }

    Collection<String> getSegmentNames() {
        segmentByName.entrySet()
    }

    String getSubSequence(String vSegment, int from, int to) {
        Util.getSubSequence(segmentByName[vSegment].sequence, from, to)
    }

    FrequencyCounter getFrequency(String vSegment) {
        countersByName[vSegment]
    }

    int getRegionSize(String vSegment, int regionId) {
        segmentByName[vSegment].regionMarkup[regionId].size()
    }

    boolean append(Clonotype clonotype) {
        def counter = countersByName[clonotype.v]
        if (counter != null) {
            counter.freq += clonotype.freq
            counter.clonotypes++
            return true
        }
        return false
    }

    @Override
    String toString() {
        "#segment\t$FrequencyCounter.HEADER\n" + segmentByName.keySet().collect {
            "${segmentByName[it]}\t${countersByName[it]}" }.join("\n")
    }
}
