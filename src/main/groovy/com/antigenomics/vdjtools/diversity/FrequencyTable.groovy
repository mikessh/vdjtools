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
 * Last modified on 2.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Countable
import com.antigenomics.vdjtools.Counter
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.pool.ClonotypeAggregator
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.sample.Sample

public class FrequencyTable {
    private final long min, max, count, diversity
    private final Map<Long, Long> frequencyMap = new HashMap<>()

    FrequencyTable(Sample sample, IntersectionType intersectionType) {
        Iterable<Countable> counters

        // collapse clonotypes by a specific key
        def clonotypeKeyGen = new ClonotypeKeyGen(intersectionType)

        def hashedCounts = new HashMap<ClonotypeKey, Counter>()

        sample.each {
            def key = clonotypeKeyGen.generateKey(it)
            def counter = hashedCounts[key]
            if (!counter)
                hashedCounts.put(key, counter = new Counter())
            counter.add(it)
        }

        this.diversity = hashedCounts.size()
        counters = hashedCounts.values()

        // compute frequency table
        long min = Long.MAX_VALUE, max = -1

        counters.each {
            long count = it.count
            frequencyMap.put(count, (frequencyMap[count] ?: 0L) + 1L)
            min = Math.min(count, min)
            max = Math.max(count, max)
        }

        this.min = min
        this.max = max
        this.count = sample.count
    }

    FrequencyTable(SampleAggregator pool) {
        long diversity = 0, count = 0

        // compute frequency table
        long min = Long.MAX_VALUE, max = -1

        pool.each { ClonotypeAggregator it ->
            long x = it.incidenceCount
            frequencyMap.put(x, (frequencyMap[x] ?: 0L) + 1L)
            min = Math.min(x, min)
            max = Math.max(x, max)
            diversity++
            count += x
        }

        this.min = min
        this.max = max
        this.count = count
        this.diversity = diversity
    }

    public long getDiversity() {
        diversity
    }

    public long getCount() {
        count
    }

    public long getAt(long clonotypeSize) {
        frequencyMap[clonotypeSize] ?: 0
    }

    public long getMin() {
        min
    }

    public long getMax() {
        max
    }

    public final List<BinInfo> getBins() {
        double s = 0
        frequencyMap.sort { it.key }.collect {
            s += it.value
            double std = Math.sqrt(it.key * (1.0 - it.key / (double) count))

            new BinInfo(it.key, Math.round(std),
                    it.key / count, std / count,
                    it.value, 1.0 - s / diversity)
        }
    }

    public static class BinInfo {
        private final long clonotypeSize, numberOfClonotypes, cloneSizeStd
        private final double complementaryCdf, clonotypeFreq, clonotypeFreqStd

        BinInfo(long clonotypeSize, long cloneSizeStd,
                double clonotypeFreq, double clonotypeFreqStd,
                long numberOfClonotypes, double complementaryCdf) {
            this.clonotypeSize = clonotypeSize
            this.cloneSizeStd = cloneSizeStd
            this.clonotypeFreq = clonotypeFreq
            this.clonotypeFreqStd = clonotypeFreqStd
            this.numberOfClonotypes = numberOfClonotypes
            this.complementaryCdf = complementaryCdf
        }

        double getCloneStd() {
            cloneSizeStd
        }

        long getClonotypeSize() {
            clonotypeSize
        }

        long getNumberOfClonotypes() {
            numberOfClonotypes
        }

        double getComplementaryCdf() {
            complementaryCdf
        }

        double getClonotypeFreq() {
            return clonotypeFreq
        }

        double getClonotypeFreqStd() {
            return clonotypeFreqStd
        }

        public static final String HEADER = "#clonotype_size\tclonotype_size_l\tclonotype_size_u\t" +
                "clonotype_freq\tclonotype_freq_l\tclonotype_freq_u\t" +
                "number_of_clonotypes\tcompl_cdf"

        @Override
        public String toString() {
            [clonotypeSize, clonotypeSize - cloneSizeStd, clonotypeSize + cloneSizeStd,
             clonotypeFreq, clonotypeFreq - clonotypeFreqStd, clonotypeFreq + clonotypeFreqStd,
             numberOfClonotypes, complementaryCdf].join("\t")
        }
    }
}
