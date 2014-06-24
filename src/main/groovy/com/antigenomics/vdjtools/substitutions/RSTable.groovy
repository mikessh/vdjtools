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

package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.Mutation
import com.antigenomics.vdjtools.segment.SegmentUtil
import com.antigenomics.vdjtools.segment.VSegmentTable

class RSTable {
    final boolean normalize
    final VSegmentTable vSegmentTable
    final Map<String, double[][][]> summaryTableBySegment = new HashMap<>()
    final double[][][] allSegmentsTable = new double[2][2][N]

    RSTable(boolean normalize, VSegmentTable vSegmentTable) {
        this.normalize = normalize
        this.vSegmentTable = vSegmentTable
        vSegmentTable.segmentNames.each {
            summaryTableBySegment.put(it, new double[2][2][N])
        }
        //summaryTableBySegment.put("All", allSegmentsTable)
    }

    void addAll(Collection<Mutation> mutations) {
        mutations.each { Mutation mutation ->
            def regionId = SegmentUtil.regionName2Id(mutation.region), vSegment = mutation.v
            def vFrequency = vSegmentTable.getFrequency(vSegment)
            def silent = mutation.isSilent ? 1 : 0

            def table = summaryTableBySegment[vSegment]
            if (table == null)
                summaryTableBySegment.put(vSegment, table = new double[2][2][N])

            def count = 1, freq = mutation.freq
            def regionSize = regionId < 5 ?
                    vSegmentTable.getRegionSize(vSegment, regionId) : mutation.parent.cdr3nt.size()

            count /= regionSize
            freq /= regionSize

            if (normalize) {
                count /= vFrequency.clonotypes
                freq /= vFrequency.freq
            }

            table[0][silent][regionId] += count
            table[1][silent][regionId] += freq
            allSegmentsTable[0][silent][regionId] += count
            allSegmentsTable[1][silent][regionId] += freq
        }
    }

    double[] rsSummary() {
        (0..<N).collect { allSegmentsTable[0][0][it] / allSegmentsTable[0][1][it] }//.collect { it.isNaN() ? 0 : it }
    }

    double[] covSummary() {
        (0..<N).collect { allSegmentsTable[0][1][it] + allSegmentsTable[0][0][it] }
    }

    double[] silentSummary() {
        (0..<N).collect { allSegmentsTable[0][1][it] }
    }

    double[] replacementSummary() {
        (0..<N).collect { allSegmentsTable[0][0][it] }
    }

    final static int N = SegmentUtil.N_REGIONS
    final static String HEADER = "Counter\t" +
            (0..<N).collect { "Frequency" }.join('\t') + '\t' +
            (0..<N).collect { "Frequency" }.join('\t') + '\t' +
            (0..<N).collect { "Clonotypes" }.join('\t') + '\t' +
            (0..<N).collect { "Clonotypes" }.join('\t') +
            "\nType\t" +
            (0..<N).collect { "Replacement" }.join('\t') + '\t' +
            (0..<N).collect { "Silent" }.join('\t') + '\t' +
            (0..<N).collect { "Replacement" }.join('\t') + '\t' +
            (0..<N).collect { "Silent" }.join('\t') +
            "\nRegion\t" +
            (1..4).collect { (0..<N).collect { SegmentUtil.regionId2Name(it) }.join("\t") }.join('\t')

    @Override
    String toString() {
        [HEADER,
         summaryTableBySegment.collect { it.key + "\t" + it.value.collect().flatten().join("\t") },
         "All\t" + allSegmentsTable.collect().flatten().join("\t")].join("\n")
    }
}
