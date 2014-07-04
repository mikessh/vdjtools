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

import com.antigenomics.vdjtools.Countable
import com.antigenomics.vdjtools.Counter
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.apache.commons.math3.util.CombinatoricsUtils

class DiversityEstimator {
    final Sample sample
    private FrequencyTable frequencyTable = null
    private DownSampler downSampler = null

    DiversityEstimator(Sample sample) {
        this.sample = sample
    }

    Diversity sampleDiversity() {
        new Diversity(sample.diversity, 0, sample.cells, false)
    }

    Diversity countNormalizedSampleDiversity(int sampleSize, int nResamples, boolean byAminoAcid) {
        downSampler = downSampler ?: new DownSampler(sample)

        if (sampleSize >= sample.cells)
            return new Diversity((long) (((double) sampleSize *
                    (double) (byAminoAcid ? sample.diversityAA : sample.diversity) / (double) sample.cells)),
                    0, sampleSize, false)

        def diversityValues = new double[nResamples]
        for (int i = 0; i < nResamples; i++) {
            def newSample = downSampler.reSample(sampleSize)
            diversityValues[i] = byAminoAcid ? newSample.diversityAA : newSample.diversity
        }

        def descrStats = new DescriptiveStatistics(diversityValues)

        return new Diversity((long) descrStats.mean, (long) descrStats.standardDeviation, sampleSize, false)
    }

    Diversity efronThisted(int maxDepth, double cvThreshold, boolean byAminoAcid) {
        this.frequencyTable = frequencyTable ?: new FrequencyTable(sample, byAminoAcid)

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
            S = sample.diversity + (double) (0..<depth).sum { int i -> h[i] * nx[i] }
            D = Math.sqrt((double) (0..<depth).sum { int i -> h[i] * h[i] * nx[i] })
            CV = D / S

            // Go to maximum count depth, but balance that STD doesn't get too high
            if (CV >= cvThreshold)
                break
        }

        new Diversity((long) S, (long) D, sample.cells, true)
    }

    Diversity chao1(boolean byAminoAcid) {
        this.frequencyTable = frequencyTable ?: new FrequencyTable(sample, byAminoAcid)
        double F1 = frequencyTable[1], F2 = frequencyTable[2], RF = F1 / F2 / 2

        new Diversity((long) (sample.diversity + F1 * RF),
                (long) (F2 * (Math.pow(RF / 2, 4) + Math.pow(2 * RF, 3) + RF * RF)),
                sample.cells, true)
    }

    private class FrequencyTable {
        private final Map<Long, Long> frequencyMap = new HashMap<>()

        FrequencyTable(Sample sample, boolean byAminoAcid) {
            Collection<Countable> counters

            if (byAminoAcid) {
                def aaCounts = new HashMap<String, Counter>()
                sample.clonotypes.each {
                    def counter = aaCounts[it.cdr3aa]
                    if (!counter)
                        aaCounts.put(it.cdr3aa, counter = new Counter())
                    counter.add(it)
                }
                counters = aaCounts.values()
            } else
                counters = sample.clonotypes

            counters.each {
                long count = it.count
                frequencyMap.put(count, (frequencyMap[count] ?: 0L) + 1L)
            }
        }

        int getAt(long count) {
            frequencyMap[count] ?: 0
        }
    }
}
