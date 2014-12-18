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

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.sample.Sample
import groovy.transform.PackageScope
import org.apache.commons.math3.util.CombinatoricsUtils

class DiversityEstimates {
    private FrequencyTable frequencyTable = null
    private final long extrapolateTo
    private final ChaoEstimator chaoEstimator

    @PackageScope
    static DiversityEstimates basicDiversityEstimates(FrequencyTable frequencyTable) {
        new DiversityEstimates(frequencyTable, -1)
    }

    public DiversityEstimates(FrequencyTable frequencyTable, long extrapolateTo) {
        this.frequencyTable = frequencyTable
        this.extrapolateTo = extrapolateTo
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    public DiversityEstimates(Sample sample, IntersectionType intersectionType, long extrapolateTo) {
        this(new FrequencyTable(sample, intersectionType), extrapolateTo)
    }

    public DiversityEstimates(SampleAggregator pool, long extrapolateTo) {
        this.frequencyTable = new FrequencyTable(pool)
        this.extrapolateTo = extrapolateTo
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    public Diversity getObservedDiversity() {
        new Diversity(frequencyTable.diversity, 0, frequencyTable.count,
                false, false, "observed_diversity")
    }

    public Diversity getEfronThisted() {
        getEfronThisted(20, 0.05)
    }

    public Diversity getEfronThisted(int maxDepth, double cvThreshold) {
        double S = -1, D = -1, CV

        for (int depth = 1; depth <= maxDepth; depth++) {
            final double[] h = new double[depth], nx = new double[depth]
            for (int y = 1; y <= depth; y++) {
                nx[y - 1] = frequencyTable[y]

                // Calculate Euler coefficients
                for (int x = 1; x <= y; x++) {
                    def coef = CombinatoricsUtils.binomialCoefficientDouble((int) (y - 1), (int) (x - 1))
                    if (x % 2 == 1)
                        h[x - 1] += coef
                    else
                        h[x - 1] -= coef
                }
            }

            // Extrapolate count
            S = frequencyTable.diversity + (double) (0..<depth).sum { int i -> h[i] * nx[i] }
            D = Math.sqrt((double) (0..<depth).sum { int i -> h[i] * h[i] * nx[i] })
            CV = D / S

            // Go to maximum count depth, but balance that STD doesn't get too high
            if (CV >= cvThreshold)
                break
        }

        new Diversity(S, D, frequencyTable.count,
                true, false, "efron_thisted")
    }

    public Diversity getChao1() {
        chaoEstimator.chao1()
    }

    public Diversity getChaoE() {
        chaoEstimator.chaoE(extrapolateTo)
    }

    public FrequencyTable getFrequencyTable() {
        frequencyTable
    }

    private static final String[] fields = ["observedDiversity", "efronThisted", "chao1", "chaoE"]
    public static final HEADER = fields.collect { [it + "_mean", it + "_std"] }.flatten().join("\t")

    @Override
    public String toString() {
        fields.collect { this."$it" }.join("\t")
    }
}
