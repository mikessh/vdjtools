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
 * A species richness estimate that reflects the total number of clonotypes in an immune repertoire.
 */
class SpeciesRichness extends DiversityEstimate {
    private final RichnessEstimateType type

    /**
     * Creates a structure holding diversity estimate summary.
     * @param mean expected value of a diversity estimate.
     * @param std standard deviation of a diversity estimate.
     * @param numberOfReads number of reads in the sample that was analyzed.
     * @param type richness estimate type. 
     */
    SpeciesRichness(long mean, long std, long numberOfReads, RichnessEstimateType type) {
        super(mean, std, numberOfReads)
        this.type = type
    }

    /**
     * Creates a structure holding diversity estimate summary.
     * @param mean expected value of a diversity estimate.
     * @param std standard deviation of a diversity estimate.
     * @param numberOfReads number of reads in the sample that was analyzed.
     * @param type richness estimate type.
     */
    SpeciesRichness(double mean, double std, long numberOfReads, RichnessEstimateType type) {
        this((long) mean, (long) std, numberOfReads, type)
    }

    /**
     * Gets the richness estimate type.
     * @return estimate type (interpolated/observed/extrapolated/lower bound estimate on total diversity).
     */
    RichnessEstimateType getType() {
        type
    }
}
