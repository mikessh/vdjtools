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
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.apache.commons.math3.util.CombinatoricsUtils
import sun.reflect.generics.reflectiveObjects.NotImplementedException

class DiversityEstimator {
    private final Sample sample
    private final IntersectionType intersectionType
    private FrequencyTable frequencyTable = null
    private DownSampler downSampler = null

    public DownSampler getDownSampler() {
        downSampler ?: (downSampler = new DownSampler(sample))
    }

    public DiversityEstimator(Sample sample, IntersectionType intersectionType) {
        this.sample = sample
        this.intersectionType = intersectionType
        this.frequencyTable = new FrequencyTable(sample, intersectionType)
    }

    public DiversityEstimator(SampleAggregator pool) {
        this.sample = null
        this.intersectionType = null
        this.frequencyTable = new FrequencyTable(pool)
    }

    public Diversity computeCollapsedSampleDiversity() {
        new Diversity(frequencyTable.diversity, 0, frequencyTable.count, false)
    }

    Diversity computeNormalizedSampleDiversity(int sampleSize, int resampleCount) {
        if (!sample)
            throw new NotImplementedException()

        if (sampleSize >= sample.count) {
            // extrapolating
            def baseDiversity = computeCollapsedSampleDiversity()
            return new Diversity((long) (((double) sampleSize * (double) baseDiversity.mean / (double) sample.count)),
                    0, sampleSize, true)
        }

        def diversityValues = new double[resampleCount]
        for (int i = 0; i < resampleCount; i++) {
            def subSample = getDownSampler().reSample(sampleSize)
            def newDiversityEstimator = new DiversityEstimator(subSample, this.intersectionType)
            diversityValues[i] = (double) newDiversityEstimator.computeCollapsedSampleDiversity().mean
        }

        def descrStats = new DescriptiveStatistics(diversityValues)

        return new Diversity((long) descrStats.mean, (long) descrStats.standardDeviation, sampleSize, false)
    }

    Diversity computeEfronThisted(int maxDepth, double cvThreshold) {
        double S = -1, D = -1, CV = -1

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

        new Diversity((long) S, (long) D, frequencyTable.count, true)
    }

    Diversity computeChao1() {
        double Sobs = frequencyTable.diversity,
               F1 = frequencyTable[1], F2 = frequencyTable[2]

        double delta = F1 * (F1 - 1) / 2 / (F2 + 1)

        new Diversity(
                (long) (Sobs + delta),
                (long) Math.sqrt(
                        delta +
                        F1 * (2 * F1 - 1) * (2 * F1 - 1) / 4 / (F2 + 1) / (F2 + 1) +
                        F1 * F1 * F2 * (F1 - 1) * (F1 - 1) / 4 / (F2 + 1) / (F2 + 1) / (F2 + 1) / (F2 + 1)
                ),
                frequencyTable.count,
                true)
    }

    FrequencyStat computeFrequencyDistributionStats() {
        new FrequencyStat(frequencyTable)
    }

    FrequencyTable getFrequencyTable() {
        frequencyTable
    }
}
