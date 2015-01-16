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

import static com.antigenomics.vdjtools.diversity.DiversityType.Observed
import static com.antigenomics.vdjtools.diversity.DiversityType.TotalDiversityLowerBoundEstimate

/**
 * Class that computes richness estimates and diversity indices. 
 * All computations are performed based on {@code FrequencyTable}.
 * NOTE: that the {@code IntersectionType} is always computed based on certain {@code IntersectionType},
 * which tells how to collapse clonotypes, e.g. consider identical CDR3 nucleotide or amino acid sequences, etc.
 * {@code IntersectionType.strict} is recommended for most purposes.
 */
public class DiversityEstimates {
    private FrequencyTable frequencyTable = null
    private final long extrapolateTo
    private final ChaoEstimator chaoEstimator

    /**
     * INTERNAL, constructor used for {@code DiversityEstimatesResampled}
     */
    @PackageScope
    static DiversityEstimates basicDiversityEstimates(FrequencyTable frequencyTable) {
        new DiversityEstimates(frequencyTable, -1)
    }

    /**
     * Creates an instance of individual-based diversity estimates class
     * @param frequencyTable a {@code FrequencyTable} summary for sample of interest
     * @param extrapolateTo extrapolated sample size. Used for ChaoE estimate. 
     *                      Typically set to the size of the largest sample when several sampels are considered.
     */
    public DiversityEstimates(FrequencyTable frequencyTable, long extrapolateTo) {
        this.frequencyTable = frequencyTable
        this.extrapolateTo = extrapolateTo
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Creates an instance of individual-based diversity estimates class. Will compute {@code FrequencyTable} for a given sample.
     * @param sample sample to be analyzed
     * @param intersectionType {@code IntersectionType} used to collapse sample during {@code FrequencyTable} computation
     * @param extrapolateTo extrapolated sample size. Used for ChaoE estimate.
     *                      Typically set to the size of the largest sample when several samples are considered.
     */
    public DiversityEstimates(Sample sample, IntersectionType intersectionType, long extrapolateTo) {
        this(new FrequencyTable(sample, intersectionType), extrapolateTo)
    }

    /**
     * Creates an instance of sample-based diversity estimates class. Will summarize clonotype occurrences to build {@code FrequencyTable}
     * @param pool a pool of several samples
     * @param extrapolateTo extrapolated sample count.
     *                      Number of samples in the largest sample pool if several are analyzed
     */
    public DiversityEstimates(SampleAggregator pool, long extrapolateTo) {
        this.frequencyTable = new FrequencyTable(pool)
        this.extrapolateTo = extrapolateTo
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Gets the observed diversity, i.e. the number of clonotypes in the sample 
     * @return
     */
    public Diversity getObservedDiversity() {
        new Diversity(frequencyTable.diversity, 0, frequencyTable.count,
                Observed, 1, "observedDiversity")
    }

    /**
     * Gets the Shannon-Weaver diversity index
     * @return
     */
    public Diversity getShannonWeaverIndex() {
        // todo: std computaiton ?
        def mean = 0, std = 0
        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            def h = bin.diversity * bin.freq * Math.log(bin.freq) // <log(p)>
            mean -= h
            //std += h * Math.log(bin.freq) // <log(p)^2>
        }
        //std -= mean * mean
        new Diversity(Math.exp(mean),
                0,//Math.sqrt((Math.exp(std) - 1) * Math.exp(std + 2 * mean)), // from lognormal approx
                frequencyTable.count,
                DiversityType.Index,
                1,
                "shannonWeaverIndex"
        )
    }

    /**
     * Gets the Inverse Simpson diversity index
     * @return
     */
    public Diversity getInverseSimpsonIndex() {
        // todo: std computaiton ?
        def mean = 0, std = 0
        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            mean += bin.diversity * bin.freq * bin.freq // <p>
            //std += bin.diversity * bin.freq * bin.freq * bin.freq // <p^2>
        }
        //std -= mean * mean
        new Diversity(1.0 / mean,
                0,
                frequencyTable.count,
                DiversityType.Index,
                1,
                "inverseSimpsonIndex"
        )
    }

    /**
     * Gets Efron-Thisted lower bound estimate of total diversity 
     * @return
     */
    public Diversity getEfronThisted() {
        getEfronThisted(20, 0.05)
    }

    /**
     * Gets Efron-Thisted lower bound estimate of total diversity
     * @param maxDepth max number of terms in series to compute
     * @param cvThreshold cv threshold, after reaching it higher order terms will be ignored
     * @return
     */
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
                TotalDiversityLowerBoundEstimate, 1, "efronThisted")
    }

    /**
     * Gets Chao1 lower bound estimate of total diversity
     * @return
     */
    public Diversity getChao1() {
        chaoEstimator.chao1()
    }

    /**
     * Gets extrapolated Chao diversity estimate 
     * @return
     */
    public Diversity getChaoE() {
        extrapolateTo < frequencyTable.count ? Diversity.DUMMY : chaoEstimator.chaoE(extrapolateTo)
    }

    /**
     * Gets the underlying frequency table
     * @return
     */
    public FrequencyTable getFrequencyTable() {
        frequencyTable
    }

    /**
     * List of fields that will be included in tabular output 
     */
    private static final String[] fields = ["observedDiversity", "chaoE",
                                            "efronThisted", "chao1",
                                            "shannonWeaverIndex", "inverseSimpsonIndex"]

    /**
     * Header string, used for tabular output
     */
    public static final HEADER = fields.collect { [it + "_mean", it + "_std"] }.flatten().join("\t")

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        fields.collect { this."$it" }.join("\t")
    }
}
