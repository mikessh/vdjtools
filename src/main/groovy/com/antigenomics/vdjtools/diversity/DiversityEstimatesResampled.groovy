/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 * Last modified on 17.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

import static com.antigenomics.vdjtools.diversity.DiversityType.Interpolated
import static com.antigenomics.vdjtools.diversity.DiversityType.TotalDiversityLowerBoundEstimate


class DiversityEstimatesResampled {
    private final Diversity observedDiversity, efronThisted, chao1
    private final int subSampleSize, resampleCount

    public DiversityEstimatesResampled(Sample sample,
                                       IntersectionType intersectionType,
                                       int subSampleSize) {
        this(sample, intersectionType, subSampleSize, 3)
    }

    public DiversityEstimatesResampled(Sample sample,
                                       IntersectionType intersectionType,
                                       int subSampleSize, int resampleCount) {
        this.subSampleSize = subSampleSize
        this.resampleCount = resampleCount

        def downSampler = new DownSampler(sample)

        def observedDiversityStat = new DescriptiveStatistics(),
            efronThistedStat = new DescriptiveStatistics(),
            chao1Stat = new DescriptiveStatistics()

        for (int i = 0; i < resampleCount; i++) {
            def subSample = downSampler.reSample(subSampleSize)
            def frequencyTable = new FrequencyTable(subSample, intersectionType)
            def diversityEstimates = DiversityEstimates.basicDiversityEstimates(frequencyTable)
            observedDiversityStat.addValue(diversityEstimates.observedDiversity.mean)
            efronThistedStat.addValue(diversityEstimates.efronThisted.mean)
            chao1Stat.addValue(diversityEstimates.chao1.mean)
        }

        this.observedDiversity = new Diversity(
                observedDiversityStat.mean,
                observedDiversityStat.standardDeviation,
                subSampleSize,
                Interpolated, true, "observed_diversity")

        this.efronThisted = new Diversity(
                efronThistedStat.mean,
                efronThistedStat.standardDeviation,
                subSampleSize,
                TotalDiversityLowerBoundEstimate, true, "efron_thisted")

        this.chao1 = new Diversity(
                chao1Stat.mean,
                chao1Stat.standardDeviation,
                subSampleSize,
                TotalDiversityLowerBoundEstimate, true, "chao1")
    }

    public Diversity getObservedDiversity() {
        observedDiversity
    }

    public Diversity getEfronThisted() {
        efronThisted
    }

    public Diversity getChao1() {
        chao1
    }

    private static final String[] fields = ["observedDiversity", "efronThisted", "chao1"]
    public static final HEADER = fields.collect { [it + "_mean", it + "_std"] }.flatten().join("\t")

    @Override
    public String toString() {
        fields.collect { this."$it" }.join("\t")
    }
}
