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
 * Last modified on 16.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.manipulation.DownSampler
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

import static com.antigenomics.vdjtools.diversity.DiversityType.*

/**
 * Class that computes richness estimates and diversity indices. 
 * Re-sampling is used to calculate the mean and standard deviation for estimates.
 * In order to normalize diversity estimates between samples, they are down-sampled to the same size,
 * typically the size of the smallest sample.
 * All computations are performed via the @see com.antigenomics.vdjtools.diversity.FrequencyTable object.
 * Note that the clonotype is always computed based on certain {@code IntersectionType},
 * which tells how to collapse clonotypes, e.g. consider identical CDR3 nucleotide or amino acid sequences, etc.
 * Therefore in some cases the result will be different from one obtained using raw clonotype frequencies.
 * This will not happen in case @see com.antigenomics.vdjtools.intersection.IntersectionType#Strict is used, 
 * which is recommended for most purposes.
 */
public class DiversityEstimatesResampled {
    private final Diversity observedDiversity, efronThisted, chao1,
                            d50Index, shannonWeaverIndex, inverseSimpsonIndex

    private final int subSampleSize, resampleCount

    /**
     * Creates an instance of individual-based diversity estimates class computed using re-sampling.
     * All computations are performed within the constructor.
     * @param sample sample to be analyzed
     * @param intersectionType {@code IntersectionType} used to collapse sample during {@code FrequencyTable} computation
     * @param subSampleSize down-sampled sample size. Typically set to the size of smallest sample if several samples are to be compared
     */
    public DiversityEstimatesResampled(Sample sample,
                                       IntersectionType intersectionType,
                                       int subSampleSize) {
        this(sample, intersectionType, subSampleSize, 3)
    }

    /**
     * Creates an instance of individual-based diversity estimates class computed using re-sampling.
     * All computations are performed within the constructor.
     * @param sample sample to be analyzed
     * @param intersectionType {@code IntersectionType} used to collapse sample during {@code FrequencyTable} computation
     * @param subSampleSize down-sampled sample size. Typically set to the size of smallest sample if several samples are to be compared
     * @param resampleCount number of re-samples to be performed
     */
    public DiversityEstimatesResampled(Sample sample,
                                       IntersectionType intersectionType,
                                       int subSampleSize, int resampleCount) {
        this.subSampleSize = subSampleSize
        this.resampleCount = resampleCount

        def downSampler = new DownSampler(sample)

        def observedDiversityStat = new DescriptiveStatistics(),
            efronThistedStat = new DescriptiveStatistics(),
            chao1Stat = new DescriptiveStatistics(),
            d50Index = new DescriptiveStatistics(),
            shannonWeaverIndexStat = new DescriptiveStatistics(),
            inverseSimpsonIndexStat = new DescriptiveStatistics()

        for (int i = 0; i < resampleCount; i++) {
            def subSample = downSampler.reSample(subSampleSize)
            def frequencyTable = new FrequencyTable(subSample, intersectionType)
            def diversityEstimates = DiversityEstimates.basicDiversityEstimates(frequencyTable)
            observedDiversityStat.addValue(diversityEstimates.observedDiversity.mean)
            efronThistedStat.addValue(diversityEstimates.efronThisted.mean)
            chao1Stat.addValue(diversityEstimates.chao1.mean)
            d50Index.addValue(diversityEstimates.d50Index.mean)
            shannonWeaverIndexStat.addValue(diversityEstimates.shannonWeaverIndex.mean)
            inverseSimpsonIndexStat.addValue(diversityEstimates.inverseSimpsonIndex.mean)
        }

        this.observedDiversity = new Diversity(
                observedDiversityStat.mean,
                observedDiversityStat.standardDeviation,
                subSampleSize,
                Interpolated, resampleCount, "observedDiversity")

        this.efronThisted = new Diversity(
                efronThistedStat.mean,
                efronThistedStat.standardDeviation,
                subSampleSize,
                TotalDiversityLowerBoundEstimate, resampleCount, "efronThisted")

        this.chao1 = new Diversity(
                chao1Stat.mean,
                chao1Stat.standardDeviation,
                subSampleSize,
                TotalDiversityLowerBoundEstimate, resampleCount, "chao1")

        this.d50Index = new Diversity(
                d50Index.mean,
                d50Index.standardDeviation,
                subSampleSize,
                Index, resampleCount, "d50index")

        this.shannonWeaverIndex = new Diversity(
                shannonWeaverIndexStat.mean,
                shannonWeaverIndexStat.standardDeviation,
                subSampleSize,
                Index, resampleCount, "shannonWeaverIndex")

        this.inverseSimpsonIndex = new Diversity(
                inverseSimpsonIndexStat.mean,
                inverseSimpsonIndexStat.standardDeviation,
                subSampleSize,
                Index, resampleCount, "inverseSimpsonIndex")
    }

    /**
     * Gets the observed diversity, i.e. the number of clonotypes in a down-sampled sample.
     * @return {@link com.antigenomics.vdjtools.diversity.Diversity} object handling the result.
     */
    public Diversity getObservedDiversity() {
        observedDiversity
    }

    /**
     * Gets Efron-Thisted lower bound estimate of total diversity computed for down-sampled clonotype frequencies.
     * @return {@link com.antigenomics.vdjtools.diversity.Diversity} object handling the result.
     */
    public Diversity getEfronThisted() {
        efronThisted
    }

    /**
     * Gets Chao1 lower bound estimate of total diversity computed for down-sampled clonotype frequencies.
     * @return {@link com.antigenomics.vdjtools.diversity.Diversity} object handling the result.
     */
    public Diversity getChao1() {
        chao1
    }

    /**
     * Gets Chao1 lower bound estimate of total diversity computed for down-sampled clonotype frequencies.
     * @return {@link com.antigenomics.vdjtools.diversity.Diversity} object handling the result.
     */
    public Diversity getD50Index() {
        d50Index
    }

    /**
     * Gets the Shannon-Weaver diversity index computed for down-sampled clonotype frequencies.
     * @return {@link com.antigenomics.vdjtools.diversity.Diversity} object handling the result.
     */
    public Diversity getShannonWeaverIndex() {
        shannonWeaverIndex
    }

    /**
     * Gets the Inverse Simpson diversity index computed for down-sampled clonotype frequencies.
     * @return {@link com.antigenomics.vdjtools.diversity.Diversity} object handling the result.
     */
    public Diversity getInverseSimpsonIndex() {
        inverseSimpsonIndex
    }

    /**
     * List of fields that will be included in tabular output 
     */
    private static final String[] fields = ["observedDiversity",
                                            "efronThisted", "chao1",
                                            "d50index", "shannonWeaverIndex", "inverseSimpsonIndex"]

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
