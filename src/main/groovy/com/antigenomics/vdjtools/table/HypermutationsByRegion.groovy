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

package com.antigenomics.vdjtools.table

import com.antigenomics.vdjtools.RegionRanges
import com.antigenomics.vdjtools.Util

class HypermutationsByRegion {
    final String species
    final Map<String, RegionRanges> regionsBySegment = new HashMap<>()
    final Map<String, double [][][]> countersBySegmentByRegion = new HashMap<>()
    final Map<String, int[]> countersBySegment = new HashMap<>()

    final static String SEP = "\\|",
                        PATTERN = /([0-9]+),(.+),([0-9]+),(.+):([ATGC])>([ATGC]),([0-9]+):(${Util.AA_LIST})>(${Util.AA_LIST})/

    HypermutationsByRegion(String species, String dataPath) {
        this.species = species
        dataPath = Util.getDataPath(dataPath)
        new File("$dataPath/regions.${species}.txt").splitEachLine('\t') {
            regionsBySegment.put(it[0], new RegionRanges(it[1], it[2],
                    it[3], it[4],
                    it[5], it[6],
                    it[7], it[8],
                    it[9], it[10]))
            countersBySegmentByRegion.put(it[0], new double[3][2][5])
            countersBySegment.put(it[0], new int[2])
        }
        countersBySegmentByRegion.put("All", new double[3][2][5])
        countersBySegment.put("All", new int[2])
    }

    void appendHypermutations(String segment, int cloneSize, String hypermutations) {
        if (!regionsBySegment.containsKey(segment)) {
            println "Unknown segment $segment for species=$species. Skipping"
        } else {
            def segmentCounters = countersBySegment[segment]
            segmentCounters[0] += cloneSize
            segmentCounters[1]++

            def segmentRegions = regionsBySegment[segment]
            def counters = countersBySegmentByRegion[segment]

            hypermutations.split(SEP).each {
                def matchList = Util.groomMatch(it =~ PATTERN)
                if (matchList) {
                    int count = matchList[0].toInteger(), pos = matchList[1].toInteger(),
                        region = segmentRegions.posToRegion(pos)
                    int silent = matchList[5] == matchList[6] ? 1 : 0

                    counters[0][silent][region] += count
                    counters[1][silent][region] += (count / cloneSize)
                    counters[2][silent][region]++
                }
            }
        }
    }

    def regions = []

    static final String HEADER = "Counter\t" +
            (0..<RegionRanges.N_REGIONS).collect { "Reads" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Reads" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Ratio" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Ratio" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Clones" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Clones" }.join('\t') +
            "\nType\t" +
            (0..<RegionRanges.N_REGIONS).collect { "Replacement" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Silent" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Replacement" }.join('\t') + '\t' +
            (0..<RegionRanges.N_REGIONS).collect { "Silent" }.join('\t') +
            "\nChain\t" +
            (0..3).collect { RegionRanges.HEADER }.join('\t')

    @Override
    String toString() {
        String table = HEADER

        countersBySegmentByRegion.each {
            table += "\n$it.key\t${it.value.collect().flatten()}"
        }

        table
    }
}
