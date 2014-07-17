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

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

class BasicStats {
    final Sample sample
    private final DescriptiveStatistics cloneSize
    private final static int SPECTRA_MIN = 8, SPECTRA_MAX = 32, SPECTRA_LEN = SPECTRA_MAX - SPECTRA_MIN + 1

    private static int spectraBin(Clonotype clonotype) {
        Math.min(SPECTRA_MAX, Math.max(SPECTRA_MIN, clonotype.cdr3aa.length())) - SPECTRA_MIN
    }
    private final double[] spectratype = new double[SPECTRA_LEN]

    BasicStats(Sample sample) {
        this.sample = sample

        cloneSize = new DescriptiveStatistics()

        sample.clonotypes.each {
            cloneSize.addValue(it.freq)
            spectratype[spectraBin(it)] += it.freq
        }
    }

    double getMeanCloneSize() {
        cloneSize.mean
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

    int getDiversityCDR3() {
        sample.diversityCDR3NT
    }

    static final String HEADER = "#cells\tdiversity\tcdr3_diversity\tmean_clone_size\tmedian_clone_size\t" +
            (SPECTRA_MIN..SPECTRA_MAX).collect { "S$it" }.join("\t")

    @Override
    String toString() {
        [cells, diversity, diversityCDR3, meanCloneSize, medianCloneSize, spectratype.collect()].flatten().join("\t")
    }
}
