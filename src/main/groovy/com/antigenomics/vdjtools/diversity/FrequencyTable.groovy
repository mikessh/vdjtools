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

import com.antigenomics.vdjtools.Countable
import com.antigenomics.vdjtools.Counter
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.sample.Sample

public class FrequencyTable {
    private final long min, max, total
    private final Map<Long, Long> frequencyMap = new HashMap<>()
    private final int diversity

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
        this.total = sample.count
    }

    public int getDiversity() {
        return diversity
    }

    public long getAt(long clonotypeSize) {
        frequencyMap[clonotypeSize] ?: 0
    }

    long getMin() {
        return min
    }

    long getMax() {
        return max
    }

    public final List<BinInfo> getBins() {
        double s = 0
        frequencyMap.sort { it.key }.collect {
            s += it.value
            new BinInfo(it.key, it.value, 1.0 - s / diversity, Math.sqrt(it.key * (1.0 - it.key / (double) total)))
        }
    }

    //@Override
    //public Iterator<BinInfo> iterator() {
    //    def innerIter = frequencyMap.entrySet().iterator()
    //    [hasNext: { innerIter.hasNext() },
    //     next   : { def entry = innerIter.next(); new BinInfo(entry.key, entry.value) }] as Iterator
    //}

    public static class BinInfo {
        private final long clonotypeSize, numberOfClonotypes
        private final double complementaryCdf, cloneSizeStd

        BinInfo(long clonotypeSize, long numberOfClonotypes, double complementaryCdf, double cloneSizeStd) {
            this.clonotypeSize = clonotypeSize
            this.numberOfClonotypes = numberOfClonotypes
            this.complementaryCdf = complementaryCdf
            this.cloneSizeStd = cloneSizeStd
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
            return complementaryCdf
        }

        public static
        final String HEADER = "#clonotype_size\tclonotype_size_l\tclonotype_size_u\tnumber_of_clonotypes\tcompl_cdf"

        @Override
        public String toString() {
            [clonotypeSize, clonotypeSize - cloneSizeStd, clonotypeSize + cloneSizeStd, numberOfClonotypes, complementaryCdf].join("\t")
        }

    }

    public static class Cdf {
        private final Map<Long, Double> cdf

        public Cdf(FrequencyTable frequencyTable) {
            def sum = 0
            frequencyTable.bins.each {

            }
        }
    }
}
