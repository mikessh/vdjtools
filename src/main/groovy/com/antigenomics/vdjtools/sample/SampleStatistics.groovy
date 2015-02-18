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
 * Last modified on 22.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.sample

/**
 * An object holding some basic statistics for a sample collection 
 */
public class SampleStatistics {
    private final long minCount, maxCount
    private final double minFreq, maxFreq
    private final int minDiversity, maxDiversity

    /**
     * Initializes a new sample statistics object with pre-computed values 
     * @param minCount minimal number of reads in a sample from given sample collection
     * @param maxCount maximal number of reads in a sample from given sample collection
     * @param minFreq minimal total frequency of clonotypes in a sample from given sample collection
     * @param maxFreq maximal total frequency of clonotypes in a sample from given sample collection
     * @param minDiversity minimal number of clonotypes in a sample from given sample collection
     * @param maxDiversity maximal number of clonotypes in a sample from given sample collection
     */
    public SampleStatistics(long minCount, long maxCount,
                            double minFreq, double maxFreq,
                            int minDiversity, int maxDiversity) {
        this.minCount = minCount
        this.maxCount = maxCount
        this.minFreq = minFreq
        this.maxFreq = maxFreq
        this.minDiversity = minDiversity
        this.maxDiversity = maxDiversity
    }

    /**
     * Gets the minimal number of reads in a sample from given sample collection
     * @return size of the smallest sample in sample collection
     */
    public long getMinCount() {
        minCount
    }

    /**
     * Gets the maximal number of reads in a sample from given sample collection
     * @return size of the largest sample in sample collection
     */
    public long getMaxCount() {
        maxCount
    }

    public double getMinFreq() {
        minFreq
    }

    public double getMaxFreq() {
        maxFreq
    }

    public int getMinDiversity() {
        minDiversity
    }

    public int getMaxDiversity() {
        maxDiversity
    }
}
