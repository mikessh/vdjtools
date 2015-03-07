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

/**
 * An estimate of repertoire diversity computed based on the abundances of sampled clonotypes.
 */
class DiversityEstimate {
    protected final mean, std
    protected final long numberOfReads

    public static DiversityEstimate DUMMY = new DiversityEstimate("NA", "NA", -1)

    /**
     * Creates a structure holding diversity estimate summary.
     * @param mean expected value of a diversity estimate.
     * @param std standard deviation of a diversity estimate.
     * @param numberOfReads number of reads in the sample that was analyzed.
     */
    DiversityEstimate(mean, std, long numberOfReads) {
        this.mean = mean
        this.std = std
        this.numberOfReads = numberOfReads
    }

    /**
     * Gets the mean value of diversity estimate. 
     * @return mean value.
     */
    def getMean() {
        mean
    }

    /**
     * Gets the standard deviation of diversity estimate. 
     * @return standard deviation.
     */
    def getStd() {
        std
    }

    /**
     * Gets the number of reads in the sample that was analyzed.  
     * @return associated read count.
     */
    long getNumberOfReads() {
        numberOfReads
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        [mean, std].join("\t")
    }
}
