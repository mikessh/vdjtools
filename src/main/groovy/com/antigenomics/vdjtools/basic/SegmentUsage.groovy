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
import com.antigenomics.vdjtools.util.ExecUtil

class SegmentUsage {
    public static boolean VERBOSE = true

    private final Map<String, double[]> vSegmentUsage = new HashMap<>(),
                                        jSegmentUsage = new HashMap<>(), vjSegmentUsage = new HashMap<>()
    private final Map<String, Integer> sampleIndex = new HashMap<>()
    private Map<String, Double> sortedVSegmTotal = new HashMap<>(), sortedJSegmTotal = new HashMap<>()
    private final int n
    private final boolean unweighted

    public SegmentUsage(SampleCollection sampleCollection, boolean unweighted) {
        this.n = sampleCollection.size()
        this.unweighted = unweighted
        // ok to use in batch intersect as all samples are pre-loaded
        sampleCollection.eachWithIndex { it, ind -> process(it, ind) }
        summarize()
    }

    public SegmentUsage(Sample[] samples, boolean unweighted) {
        this.n = samples.length
        this.unweighted = unweighted
        // ok to use in batch intersect as all samples are pre-loaded
        samples.eachWithIndex { it, ind -> process(it, ind) }
        summarize()
    }

    private void process(Sample sample, int index) {
        ExecUtil.report(this, "[${new Date()} SegmentUsage] Processing sample ${sample.sampleMetadata.sampleId}", VERBOSE)
        sample.each { Clonotype clonotype ->
            def vArray = vSegmentUsage[clonotype.v],
                jArray = jSegmentUsage[clonotype.j],
                vjArray = vjSegmentUsage[clonotype.v + "\t" + clonotype.j]

            if (!vArray)
                vSegmentUsage.put(clonotype.v, vArray = new double[n])

            if (!jArray)
                jSegmentUsage.put(clonotype.j, jArray = new double[n])

            if (!vjArray)
                vjSegmentUsage.put(clonotype.v + "\t" + clonotype.j, vjArray = new double[n])

            def increment = unweighted ? 1 : clonotype.freq

            vArray[index] += increment
            jArray[index] += increment
            vjArray[index] += increment
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
        jUsageVector(getSampleIndex(sampleId))
    }

    public double[] jUsageVector(int sampleIndex) {
        usageVector(jSegmentUsage, sortedJSegmTotal, sampleIndex)
    }

    public double[] vUsageVector(String sampleId) {
        vUsageVector(getSampleIndex(sampleId))
    }

    public double[] vUsageVector(int sampleIndex) {
        usageVector(vSegmentUsage, sortedVSegmTotal, sampleIndex)
    }

    private int getSampleIndex(String sampleId) {
        if (!this.sampleIndex.containsKey(sampleId))
            throw new IllegalArgumentException("$sampleId is not in the sample collection used to build segment usage matrix")
        this.sampleIndex[sampleId]
    }

    private static double[] usageVector(Map<String, double[]> usageMap, Map<String, Double> totalMap, int sampleIndex) {
        def sampleTotal = usageMap.values().collect { it[sampleIndex] }.sum()
        totalMap.collect {
            usageMap[it.key][sampleIndex] / (double) (sampleTotal + 1e-7)
        } as double[]
    }

    public double[][] vjUsageMatrix(String sampleId) {
        vjUsageMatrix(getSampleIndex(sampleId))
    }

    public double[][] vjUsageMatrix(int sampleIndex) {
        double sampleTotal = (double) jSegmentUsage.values().collect { it[sampleIndex] }.sum()

        def matrix = new double[jSegmentUsage.size()][vSegmentUsage.size()]

        sortedJSegmTotal.eachWithIndex { jEntry, ii ->
            sortedVSegmTotal.eachWithIndex { vEntry, jj ->
                def v = vEntry.key, j = jEntry.key, vj = v + "\t" + j
                matrix[ii][jj] = vjSegmentUsage.containsKey(vj) ? (vjSegmentUsage[vj][sampleIndex] / sampleTotal) : 0d
            }
        }

        matrix
    }

    public String[] jUsageHeader() {
        sortedJSegmTotal.collect { it.key }
    }

    public String[] vUsageHeader() {
        sortedVSegmTotal.collect { it.key }
    }
}
