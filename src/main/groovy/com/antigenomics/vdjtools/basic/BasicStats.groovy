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

import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

class BasicStats {
    final Sample sample
    private final DescriptiveStatistics cloneSize, cdr3ntLength, insertSize, ndnSize
    private int oofCount
    private double oofRatio

    BasicStats(Sample sample) {
        this.sample = sample

        this.cloneSize = new DescriptiveStatistics()
        this.cdr3ntLength = new DescriptiveStatistics()
        this.insertSize = new DescriptiveStatistics()
        this.ndnSize = new DescriptiveStatistics()

        sample.each {
            cloneSize.addValue(it.freq)
            cdr3ntLength.addValue(it.cdr3nt.length())
            insertSize.addValue(it.insertSize)
            ndnSize.addValue(it.NDNSize)
            if (!it.inFrame) {
                oofCount++
                oofRatio += it.freq
            }
        }
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

    long getCells() {
        sample.count
    }

    int getDiversity() {
        sample.diversity
    }

    static final String HEADER = "cells\tdiversity\t" +
            "mean_clone_fraction\tmedian_clone_fraction\t" +
            "oof_count\toof_fraction\t" +
            "mean_cdr3nt_length\tmean_insert_size\tmean_ndn_size"

    @Override
    String toString() {
        [cells, diversity,
         meanCloneSize, medianCloneSize,
         oofCount, oofRatio,
         meanCdr3ntLength, meanInsertSize, meanNDNSize
        ].flatten().join("\t")
    }
}
