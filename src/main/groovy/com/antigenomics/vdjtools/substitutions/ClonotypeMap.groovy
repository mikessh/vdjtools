package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Mutation
import com.antigenomics.vdjtools.igblast.ClonotypeParsedData
import com.antigenomics.vdjtools.igblast.MutationParseData
import com.antigenomics.vdjtools.segment.VSegmentTable

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
class ClonotypeMap {
    final Map<String, ClonotypeParsedData> innerMap = new HashMap()
    final Map<String, List<ClonotypeParsedData>> byCdr3Map = new HashMap()
    final Map<String, Map<MutationParseData, AlleleCounter>> freqByMutByV = new HashMap()
    final VSegmentTable vSegmentTable

    public ClonotypeMap(String species, String fileName) {
        vSegmentTable = new VSegmentTable(species)

        new File(fileName).eachLine { line ->
            if (!line.startsWith("#")) {
                def clonotype = new ClonotypeParsedData(line)

                if (innerMap.containsKey(clonotype.key)) {
                    println "[WARNING] Duplicate clonotype (identical CDR3 + mutations) found: " +
                            "${clonotype.displayName}. Further considering only the one with highest count."

                    innerMap[clonotype.key].count += clonotype.count
                    innerMap[clonotype.key].freq += clonotype.freq
                } else {
                    boolean goodV = vSegmentTable.append(clonotype)

                    if (goodV) {
                        innerMap.put(clonotype.key, clonotype)

                        def freqByMut = freqByMutByV[clonotype.v]

                        if (freqByMut == null)
                            freqByMutByV.put(clonotype.v, freqByMut = new HashMap<MutationParseData, AlleleCounter>())

                        clonotype.mutations.each { MutationParseData mpd ->
                            def mc = freqByMut[mpd]
                            if (mc == null)
                                freqByMut.put(mpd, mc = new AlleleCounter())
                            mc.freq = mc.freq + clonotype.freq
                            mc.cdr3Len.add(clonotype.cdr3nt.length())
                        }

                        def clonotypeList = byCdr3Map[clonotype.cdr3nt]
                        if (clonotypeList == null)
                            byCdr3Map.put(clonotype.cdr3nt, clonotypeList = new ArrayList<ClonotypeParsedData>())
                        clonotypeList.add(clonotype)
                    } else {
                        println "[WARNING] Unrecognized V found: ${clonotype.displayName}. Skipping"
                    }
                }
            }
        }
    }

    Collection<ClonotypeParsedData> getClonotypes() {
        innerMap.values()
    }

    void deduceAlleles(double freqThreshold, double vFreqThreshold, int spectraThreshold) {
        clonotypes.each { clonotype ->
            def vFreq = vSegmentTable.getFrequency(clonotype).freq,
                freqByMut = freqByMutByV[clonotype.v]

            clonotype.mutations.each { mpd ->
                def mutCounter = freqByMut[mpd]

                if (vFreq < vFreqThreshold ||
                        mutCounter.freq / vFreq < freqThreshold ||
                        mutCounter.cdr3Len.size() < spectraThreshold)
                    clonotype.shms.add(mpd)
                else
                    clonotype.alleles.add(mpd)
            }
        }
    }

    private class AlleleCounter {
        double freq = 0
        final Set<Integer> cdr3Len = new HashSet<>()
    }
}
