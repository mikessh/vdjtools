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
 * Last modified on 4.3.2015 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.sample.Sample
import groovy.transform.PackageScope
import org.apache.commons.math3.util.CombinatoricsUtils

/**
 * Class that computes richness estimates and diversity indices. 
 * All computations are performed via the {@link com.antigenomics.vdjtools.diversity.FrequencyTable} object.
 * Note that the clonotype is always computed based on certain {@link com.antigenomics.vdjtools.intersection.IntersectionType},
 * which tells how to collapse clonotypes, e.g. consider identical CDR3 nucleotide or amino acid sequences, etc.
 * Therefore in some cases the result will be different from one obtained using raw clonotype frequencies.
 * This will not happen in case {@link com.antigenomics.vdjtools.intersection.IntersectionType#Strict} is used,
 * which is recommended for most purposes.
 *
 * @see ResamplingEstimator
 * @see com.antigenomics.vdjtools.diversity.DiversityEstimate
 */
class ExactEstimator extends DiversityEstimator {
    private final long extrapolateTo
    private final ChaoEstimator chaoEstimator

    /**
     * INTERNAL, constructor used for {@link ResamplingEstimator.}
     */
    @PackageScope
    static ExactEstimator basicDiversityEstimates(FrequencyTable frequencyTable) {
        new ExactEstimator(frequencyTable, -1)
    }

    /**
     * Creates an instance of individual-based diversity estimates class.
     * @param frequencyTable a {@link com.antigenomics.vdjtools.diversity.FrequencyTable} summary for sample of interest.
     * @param extrapolateTo desired extrapolation extent. Used to compute the
     * {@link com.antigenomics.vdjtools.diversity.ChaoEstimator#chaoE} estimate.
     *                      For most cases, it should be set to the size of the largest sample when several samples are to be compared.
     */
    ExactEstimator(FrequencyTable frequencyTable, long extrapolateTo) {
        super(frequencyTable, EstimationMethod.Exact)
        this.extrapolateTo = extrapolateTo
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Creates an instance of individual-based diversity estimates class.
     * Will compute {@link com.antigenomics.vdjtools.diversity.FrequencyTable} for a given sample.
     * @param sample a sample to analyze.
     * @param intersectionType {@link com.antigenomics.vdjtools.intersection.IntersectionType} that will be used
     *                         to collapse sample during {@link com.antigenomics.vdjtools.diversity.FrequencyTable} computation
     * @param extrapolateTo desired extrapolation extent. Used to compute the
     * {@link com.antigenomics.vdjtools.diversity.ChaoEstimator#chaoE} estimate.
     *                      For most cases, it should be set to the size of the largest sample when several samples are to be compared.
     */
    ExactEstimator(Sample sample, IntersectionType intersectionType, long extrapolateTo) {
        this(new FrequencyTable(sample, intersectionType), extrapolateTo)
    }

    /**
     * Creates an instance of sample-based diversity estimates class.
     * Will summarize clonotype occurrences to build {@link com.antigenomics.vdjtools.diversity.FrequencyTable}
     * @param pool a pool of several samples to analyze.
     * @param extrapolateTo desired extrapolation extent. Used to compute the
     * {@link com.antigenomics.vdjtools.diversity.ChaoEstimator#chaoE} estimate.
     *                      For most cases, it should be set to the number of samples in the largest sample pool if several are to be compared.
     */
    ExactEstimator(SampleAggregator pool, long extrapolateTo) {
        this(new FrequencyTable(pool), extrapolateTo)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getShannonWeinerIndex() {
        // todo: std computaiton ?
        def mean = 0
        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            def h = bin.diversity * bin.freq * Math.log(bin.freq)
            mean -= h
        }
        new DiversityIndex(Math.exp(mean), 0, frequencyTable.count)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getInverseSimpsonIndex() {
        // todo: std computaiton ?
        def mean = 0
        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            mean += bin.diversity * bin.freq * bin.freq
        }
        //std -= mean * mean
        new DiversityIndex(1.0 / mean, 0, frequencyTable.count)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getDxxIndex(double fraction) {
        if (fraction < 0 || fraction > 1)
            throw new IllegalArgumentException("Fraction value should be within [0,1] bounds.")

        def div = 0, freqSum = 0

        frequencyTable.bins.find {
            freqSum += it.freq
            div++

            freqSum >= fraction
        }

        new DiversityIndex(div / (double) frequencyTable.diversity, 0, frequencyTable.count)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getEfronThisted(int maxDepth, double cvThreshold) {
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

        new SpeciesRichness((long) S, (long) D, frequencyTable.count, RichnessEstimateType.TotalDiversityLowerBoundEstimate)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getChao1() {
        chaoEstimator.chao1()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityEstimate getChaoE() {
        extrapolateTo < frequencyTable.count ? DiversityEstimate.DUMMY : chaoEstimator.chaoE(extrapolateTo)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getObservedDiversity() {
        new SpeciesRichness(frequencyTable.diversity, 0, frequencyTable.count, RichnessEstimateType.Observed)
    }

}
