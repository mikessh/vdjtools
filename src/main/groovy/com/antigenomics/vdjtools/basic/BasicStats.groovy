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
    private int oofCount
    private double oofRatio
    private final long count
    private final int diversity
    private final double convergence

    BasicStats(Sample sample) {
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
            if (!it.inFrame) {
                oofCount++
                oofRatio += it.freq
            }

            aaSet.add(clonotypeKeyGen.generateKey(it))
        }

        this.convergence = diversity / (double) aaSet.size() - 1
    }

    double getMeanCloneSize() {
        cloneSize.mean
    }

    double getMeanCdr3ntLength() {
        cdr3ntLength.mean
    }

    double getMeanNDNSize() {
        ndnSize.mean
    }

    double getMeanInsertSize() {
        insertSize.mean
    }

    double getMedianCloneSize() {
        cloneSize.getPercentile(50)
    }

    long getCount() {
        count
    }

    int getDiversity() {
        diversity
    }

    double getConvergence() {
        return convergence
    }

    static final String HEADER = "count\tdiversity\t" +
            "mean_clone_fraction\tmedian_clone_fraction\t" +
            "oof_count\toof_fraction\t" +
            "mean_cdr3nt_length\tmean_insert_size\tmean_ndn_size\t" +
            "convergence"

    @Override
    String toString() {
        [count, diversity,
         meanCloneSize, medianCloneSize,
         oofCount, oofRatio,
         meanCdr3ntLength, meanInsertSize, meanNDNSize,
         convergence
        ].flatten().join("\t")
    }
}
