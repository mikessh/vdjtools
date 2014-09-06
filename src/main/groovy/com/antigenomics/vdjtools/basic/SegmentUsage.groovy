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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection

class SegmentUsage {
    private final Map<String, double[]> vSegmentUsage = new HashMap<>(), jSegmentUsage = new HashMap<>()
    private final Map<String, Integer> sampleIndex = new HashMap<>()
    private Map<String, Double> sortedVSegmTotal = new HashMap<>(), sortedJSegmTotal = new HashMap<>()
    private final sampleCollection
    private final boolean unweighted

    public SegmentUsage(SampleCollection sampleCollection, boolean unweighted) {
        this.sampleCollection = sampleCollection
        this.unweighted = unweighted
        // ok to use in batch intersect as all samples are pre-loaded
        sampleCollection.eachWithIndex { it, ind -> process(it, ind) }
        summarize()
    }

    private void process(Sample sample, int index) {
        println "[${new Date()} SegmentUsage] Processing sample ${sample.sampleMetadata.sampleId}"
        sample.each { Clonotype clonotype ->
            def vArray = vSegmentUsage[clonotype.v],
                jArray = jSegmentUsage[clonotype.j]

            if (!vArray)
                vSegmentUsage.put(clonotype.v, vArray = new double[sampleCollection.size()])

            if (!jArray)
                jSegmentUsage.put(clonotype.j, jArray = new double[sampleCollection.size()])

            def increment = unweighted ? 1 : clonotype.freq

            vArray[index] += increment
            jArray[index] += increment
        }
        sampleIndex.put(sample.sampleMetadata.sampleId, index)
    }

    private void summarize() {
        vSegmentUsage.each {
            sortedVSegmTotal.put(it.key, (double) it.value.collect().sum() ?: 0)
        }
        sortedVSegmTotal = sortedVSegmTotal.sort()
        jSegmentUsage.each {
            sortedJSegmTotal.put(it.key, (double) it.value.collect().sum() ?: 0)
        }
        sortedJSegmTotal = sortedJSegmTotal.sort()
    }

    public double[] jUsageVector(String sampleId) {
        if (!sampleIndex.containsKey(sampleId))
            throw new IllegalArgumentException("$sampleId is not in the sample collection used to build usage matrix")
        def index = sampleIndex[sampleId]
        def sampleTotal = jSegmentUsage.values().collect { it[index] }.sum()
        sortedJSegmTotal.collect {
            jSegmentUsage[it.key][index] / (sampleTotal + 1e-7)
        } as double[]
    }

    public double[] vUsageVector(String sampleId) {
        if (!sampleIndex.containsKey(sampleId))
            throw new IllegalArgumentException("$sampleId is not in the sample collection used to build usage matrix")
        def index = sampleIndex[sampleId]
        def sampleTotal = jSegmentUsage.values().collect { it[index] }.sum()
        sortedVSegmTotal.collect {
            vSegmentUsage[it.key][index] / (sampleTotal + 1e-7)
        } as double[]
    }

    public double vJSD(String sampleId1, String sampleId2) {
        double[] p = vUsageVector(sampleId1), q = vUsageVector(sampleId2)

        (0..<p.length).collect { int i ->
            double m = (p[i] + q[i]) / 2.0
            (p[i] > 0 ? (Math.log(p[i] / m) * p[i]) : 0d) + (q[i] > 0 ? (Math.log(q[i] / m) * q[i]) : 0d)
        }.sum() / 2.0 / Math.log(2.0)
    }

    public String[] jUsageHeader() {
        sortedJSegmTotal.collect { it.key }
    }

    public String[] vUsageHeader() {
        sortedVSegmTotal.collect { it.key }
    }
}
