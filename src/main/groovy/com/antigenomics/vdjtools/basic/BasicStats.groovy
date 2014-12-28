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

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

class BasicStats {
    private final DescriptiveStatistics cloneSize, cdr3ntLength, insertSize, ndnSize
    private int ncDiversity
    private double ncFrequency
    private final long count
    private final int diversity
    private final double convergence

    public BasicStats(Sample sample) {
        this.count = sample.count
        this.diversity = sample.diversity

        this.cloneSize = new DescriptiveStatistics()
        this.cdr3ntLength = new DescriptiveStatistics()
        this.insertSize = new DescriptiveStatistics()
        this.ndnSize = new DescriptiveStatistics()

        def clonotypeKeyGen = new ClonotypeKeyGen(IntersectionType.AminoAcidVJ)
        def aaSet = new HashSet<ClonotypeKey>()

        sample.each {
            cloneSize.addValue(it.freq)
            cdr3ntLength.addValue(it.cdr3nt.length())

            def x = it.insertSize
            if (x > -1)
                insertSize.addValue(x)

            ndnSize.addValue(it.NDNSize)
            if (!it.coding) {
                ncDiversity++
                ncFrequency += it.freq
            }

            aaSet.add(clonotypeKeyGen.generateKey(it))
        }

        this.convergence = diversity / (double) aaSet.size() - 1
    }

    public double getMeanFrequency() {
        cloneSize.mean
    }

    public double getMeanCdr3ntLength() {
        cdr3ntLength.mean
    }

    public double getMeanNDNSize() {
        ndnSize.mean
    }

    public double getMeanInsertSize() {
        insertSize.mean
    }

    public double getMedianFrequency() {
        cloneSize.getPercentile(50)
    }

    public long getCount() {
        count
    }

    public int getDiversity() {
        diversity
    }

    public int getNcDiversity() {
        ncDiversity
    }

    public double getNcFrequency() {
        ncFrequency
    }

    public double getConvergence() {
        return convergence
    }

    static final String HEADER = "count\tdiversity\t" +
            "mean_frequency\tmedian_frequency\t" +
            "nc_diversity\tnc_frequency\t" +
            "mean_cdr3nt_length\tmean_insert_size\tmean_ndn_size\t" +
            "convergence"

    @Override
    String toString() {
        [count, diversity,
         meanFrequency, medianFrequency,
         ncDiversity, ncFrequency,
         meanCdr3ntLength, meanInsertSize, meanNDNSize,
         convergence
        ].flatten().join("\t")
    }
}
