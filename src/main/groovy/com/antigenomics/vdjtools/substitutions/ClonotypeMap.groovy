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
 *
 * Last modified on 19.6.2014 by mikesh
 */

package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Mutation
import com.antigenomics.vdjtools.segment.VSegmentTable

class ClonotypeMap {
    private final Map<String, Clonotype> innerMap = new HashMap()
    private final Map<String, List<Clonotype>> byCdr3Map = new HashMap()
    private final Map<String, Map<Mutation, AlleleCounter>> freqByMutByV = new HashMap()
    private final VSegmentTable vSegmentTable

    public ClonotypeMap(VSegmentTable vSegmentTable, Collection<Clonotype> clonotypes) {
        this.vSegmentTable = vSegmentTable

        clonotypes.each { clonotype ->
            if (innerMap.containsKey(clonotype.key)) {
                println "[WARNING] Duplicate clonotype (identical CDR3 + mutations) found: " +
                        "${clonotype.displayName}. Appending count and skipping."

                innerMap[clonotype.key].count += clonotype.count
                innerMap[clonotype.key].freq += clonotype.freq
            } else {
                boolean goodV = vSegmentTable.append(clonotype)

                if (goodV) {
                    innerMap.put(clonotype.key, clonotype)

                    def freqByMut = freqByMutByV[clonotype.v]

                    if (freqByMut == null)
                        freqByMutByV.put(clonotype.v, freqByMut = new HashMap<Mutation, AlleleCounter>())

                    clonotype.mutations.each { Mutation mpd ->
                        def mc = freqByMut[mpd]
                        if (mc == null)
                            freqByMut.put(mpd, mc = new AlleleCounter())
                        mc.freq += clonotype.freq
                        mc.cdr3Len.add(clonotype.cdr3nt.length())
                    }

                    def clonotypeList = byCdr3Map[clonotype.cdr3nt]
                    if (clonotypeList == null)
                        byCdr3Map.put(clonotype.cdr3nt, clonotypeList = new ArrayList<Clonotype>())
                    clonotypeList.add(clonotype)
                } else {
                    println "[WARNING] Unrecognized V $clonotype.v found in $clonotype.displayName. Skipping"
                }
            }
        }
    }

    Collection<Clonotype> getClonotypes() {
        innerMap.values()
    }

    Collection<List<Clonotype>> getClonotypesByCdr3() {
        byCdr3Map.values()
    }

    List<Clonotype> getByCdr3(String cdr3nt) {
        byCdr3Map[cdr3nt]
    }

    Clonotype getByKey(String key) {
        innerMap[key]
    }

    void deduceAlleles(double freqThreshold, int spectraThreshold) {
        clonotypes.each { clonotype ->
            String vSegment = clonotype.v
            def vFreq = vSegmentTable.getFrequency(vSegment).freq,
                freqByMut = freqByMutByV[vSegment]

            clonotype.mutations.each { mpd ->
                def mutCounter = freqByMut[mpd]

                if (mutCounter.freq / vFreq < freqThreshold ||
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
