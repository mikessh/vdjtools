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
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.manipulation.DownSampler
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import sun.reflect.generics.reflectiveObjects.NotImplementedException

/**
 * Class that computes richness estimates and diversity indices. 
 * Re-sampling is used to calculate the mean and standard deviation for estimates.
 * In order to normalize diversity estimates between samples, they are down-sampled to the same size,
 * typically the size of the smallest sample.
 * All computations are performed via the {@link com.antigenomics.vdjtools.diversity.FrequencyTable} object.
 * Note that the clonotype is always computed based on certain {@link OverlapType},
 * which tells how to collapse clonotypes, e.g. consider identical CDR3 nucleotide or amino acid sequences, etc.
 * Therefore in some cases the result will be different from one obtained using raw clonotype frequencies.
 * This will not happen in case {@link OverlapType#Strict} is used,
 * which is recommended for most purposes.
 */
class ResamplingEstimator extends DiversityEstimator {
    private final DiversityIndex d50Index, shannonWeinerIndex, inverseSimpsonIndex
    private final SpeciesRichness observedDiversity, efronThisted, chao1

    protected final int subSampleSize, resampleCount

    /**
     * Creates an instance of individual-based diversity estimates class computed using re-sampling.
     * All computations are performed within the constructor.
     * @param sample sample to be analyzed
     * @param intersectionType {@code IntersectionType} used to collapse sample during {@code FrequencyTable} computation
     * @param subSampleSize down-sampled sample size. Typically set to the size of smallest sample if several samples are to be compared
     */
    ResamplingEstimator(Sample sample,
                        OverlapType intersectionType,
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
    ResamplingEstimator(Sample sample,
                        OverlapType intersectionType,
                        int subSampleSize, int resampleCount) {
        super(null, EstimationMethod.Resampled)
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
            def diversityEstimates = ExactEstimator.basicDiversityEstimates(frequencyTable)
            observedDiversityStat.addValue(diversityEstimates.observedDiversity.mean)
            efronThistedStat.addValue(diversityEstimates.efronThisted.mean)
            chao1Stat.addValue(diversityEstimates.chao1.mean)
            d50Index.addValue(diversityEstimates.d50Index.mean)
            shannonWeaverIndexStat.addValue(diversityEstimates.shannonWeinerIndex.mean)
            inverseSimpsonIndexStat.addValue(diversityEstimates.inverseSimpsonIndex.mean)
        }

        this.d50Index = new DiversityIndex(
                d50Index.mean,
                d50Index.standardDeviation,
                subSampleSize)

        this.shannonWeinerIndex = new DiversityIndex(
                shannonWeaverIndexStat.mean,
                shannonWeaverIndexStat.standardDeviation,
                subSampleSize)

        this.inverseSimpsonIndex = new DiversityIndex(
                inverseSimpsonIndexStat.mean,
                inverseSimpsonIndexStat.standardDeviation,
                subSampleSize)

        this.efronThisted = new SpeciesRichness(
                efronThistedStat.mean,
                efronThistedStat.standardDeviation,
                subSampleSize,
                RichnessEstimateType.TotalDiversityLowerBoundEstimate)

        this.chao1 = new SpeciesRichness(
                chao1Stat.mean,
                chao1Stat.standardDeviation,
                subSampleSize,
                RichnessEstimateType.TotalDiversityLowerBoundEstimate)

        this.observedDiversity = new SpeciesRichness(
                observedDiversityStat.mean,
                observedDiversityStat.standardDeviation,
                subSampleSize,
                RichnessEstimateType.Observed)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getDxxIndex(double fraction) {
        throw new NotImplementedException()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getD50Index() {
        d50Index
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getShannonWeinerIndex() {
        shannonWeinerIndex
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityIndex getInverseSimpsonIndex() {
        inverseSimpsonIndex
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getObservedDiversity() {
        observedDiversity
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getEfronThisted(int maxDepth, double cvThreshold) {
        throw new NotImplementedException()
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getEfronThisted() {
        efronThisted
    }

    /**
     * {@inheritDoc}
     */
    @Override
    SpeciesRichness getChao1() {
        chao1
    }

    /**
     * {@inheritDoc}
     */
    @Override
    DiversityEstimate getChaoE() {
        DiversityEstimate.DUMMY
    }

    @Override
    FrequencyTable getFrequencyTable() {
        throw new NotImplementedException()
    }
}
