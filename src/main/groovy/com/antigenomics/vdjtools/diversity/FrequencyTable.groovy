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

/**
 * A base class for providing info on the frequencies of rare (singletons, doubletons)
 * and abundant clonotypes.
 * The class is a wrapper for frequency -> number of clonotypes with a given frequency table
 */
public class FrequencyTable {
    private final long count, diversity
    private final Map<Long, FrequencyTableBin> frequencyMap = new HashMap<>()

    /**
     * Creates a frequency table from cache
     * @param frequencyTableCache
     */
    public FrequencyTable(Map<Long, Long> frequencyTableCache) {
        long count = 0, diversity = 0
        frequencyTableCache.each {
            frequencyMap.put(it.key, new FrequencyTableBin(it.key, it.value))
            count += it.key * it.value
            diversity += it.value
        }
        this.count = count
        this.diversity = diversity
    }

    /**
     * Gets the frequency table cache
     * @return
     */
    public Map<Long, Long> getCache() {
        frequencyMap.values().collectEntries {
            [(it.count): it.diversity]
        }
    }

    /**
     * Creates frequency table that bins clonotypes according to their frequency
     * @param sample sample to bin
     * @param intersectionType intersection type used to collapse clonotypes
     */
    public FrequencyTable(Sample sample, IntersectionType intersectionType) {
        Iterable<Countable> counters

        // collapse clonotypes by a specific key
        def clonotypeKeyGen = new ClonotypeKeyGen(intersectionType)

        // todo: parallelization possible
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
        counters.each {
            long x = it.count
            def bin = frequencyMap[x]
            if (!bin)
                frequencyMap.put(x, bin = new FrequencyTableBin(x))
            bin.increment()
        }

        this.count = sample.count
    }

    /**
     * Creates frequency table that bins clonotypes according to their occurrence in a set of samples
     * @param pool pooled samples
     */
    FrequencyTable(SampleAggregator pool) {
        long diversity = 0, count = 0

        // compute frequency table
        pool.each { ClonotypeAggregator it ->
            long x = it.incidenceCount
            def bin = frequencyMap[x]
            if (!bin)
                frequencyMap.put(x, bin = new FrequencyTableBin(x))
            bin.increment()
            diversity++
            count += x
        }

        this.count = count
        this.diversity = diversity
    }

    /**
     * Gets the total number of clonotypes in this table
     * (in accordance with IntersectionType used to collapse the sample)
     * @return
     */
    public long getDiversity() {
        diversity
    }

    /**
     * Gets the total read count in this table
     * @return
     */
    public long getCount() {
        this.count
    }

    /**
     * Gets the number of singleton clonotypes, i.e. those met only once
     * @return
     */
    public int getSingletons() {
        this[1]
    }

    /**
     * Gets the number of doubleton clonotypes, i.e. those met only twice
     * @return
     */
    public int getDoubletons() {
        this[2]
    }

    /**
     * Gets the number of clonotypes with the given count
     * @param count number of reads
     * @return
     */
    public int getAt(long count) {
        def bin = frequencyMap[count]
        bin ? bin.diversity : 0
    }

    /**
     * Gets the bins that have at least one clonotype in them
     * @return
     */
    public Collection<FrequencyTableBin> getBins() {
        Collections.unmodifiableCollection(frequencyMap.values())
    }

    /**
     * Frequency table bin
     */
    public class FrequencyTableBin {
        private final long count
        private int diversity = 0

        /**
         * Creates a new frequency table bin
         * @param count number of reads, used as key
         */
        public FrequencyTableBin(long count) {
            this.count = count
        }

        /**
         * Creates a new frequency table bin
         * @param count number of reads, used as key
         * @param diversity number of clonotypes in this bin
         */
        public FrequencyTableBin(long count, long diversity) {
            this.count = count
            this.diversity = diversity
        }

        /**
         * Increment the number of clonotypes in this bin
         */
        public void increment() {
            this.diversity++
        }

        /**
         * Decrement the number of clonotypes in this bin
         */
        public void decrement() {
            this.diversity--
        }

        /**
         * Get the number of reads that specify this bin
         * @return
         */
        public long getCount() {
            this.count
        }

        /**
         * Gets the number of clonotypes in this bin
         * @return
         */
        public int getDiversity() {
            this.diversity
        }

        /**
         * Gets the ratio of clonotypes in this bin and the total diversity of the sample
         * (in accordance with IntersectionType used to collapse the sample)
         * @return
         */
        public double getRelativeDiversity() {
            this.diversity / (double) FrequencyTable.this.diversity
        }

        /**
         * Get the clonotype frequency that specify this bin
         * @return
         */
        public double getFreq() {
            this.count / (double) FrequencyTable.this.count
        }

        @Override
        boolean equals(o) {
            if (this.is(o)) return true
            if (getClass() != o.class) return false

            FrequencyTableBin that = (FrequencyTableBin) o

            this.count == that.count
        }

        @Override
        int hashCode() {
            (int) (this.count ^ (this.count >>> 32))
        }
    }
}
