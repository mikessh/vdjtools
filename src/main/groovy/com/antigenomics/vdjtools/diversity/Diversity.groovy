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

/**
 * Class holding general information on a diversity estimate
 */
class Diversity {
    public final long mean, std, n
    public final int resamples
    public final DiversityType type
    public final String estimatorName

    /**
     * Dummy value used to represent diversity estimate that could not be computed,
     * e.g. counting ChaoE for {@code n} smaller than sample size
     */
    public static Diversity DUMMY = new Diversity(-1, -1, -1, DiversityType.Unknown, -1, "DUMMY")

    /**
     * Creates a structure holding diversity estimate summary. Double precision constructor
     * @param mean expected value of a diversity estimate
     * @param std standard deviaiton of a diversity estimate
     * @param n number of reads in the sample that was analyzed
     * @param type diversity estimate type
     * @param resamples number of re-samples used to compute this estimate
     * @param estimatorName estimator short name
     */
    Diversity(double mean, double std, long n,
              DiversityType type, int resamples, String estimatorName) {
        this((long) mean, (long) std, n, type, resamples, estimatorName)
    }

    /**
     * Creates a structure holding diversity estimate summary. Int64 precision constructor.
     * @param mean expected value of a diversity estimate
     * @param std standard deviaiton of a diversity estimate
     * @param n number of reads in the sample that was analyzed
     * @param type diversity estimate type
     * @param resamples number of re-samples used to compute this estimate
     * @param estimatorName estimator short name
     */
    Diversity(long mean, long std, long n,
              DiversityType type, int resamples, String estimatorName) {
        this.mean = mean
        this.std = std
        this.n = n
        this.resamples = resamples
        this.type = type
        this.estimatorName = estimatorName
        this.HEADER = [estimatorName + "_mean", estimatorName + "_std"]
    }

    /**
     * Header string, used for tabular output
     */
    public final String HEADER
    
    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        [mean, std].join("\t")
    }
}