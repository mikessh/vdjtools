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
 * Type of richness estimate, mainly used by {@link com.antigenomics.vdjtools.diversity.Rarefaction}
 */
enum RichnessEstimateType {
    Interpolated(0), Observed(1), Extrapolated(2), TotalDiversityLowerBoundEstimate(3)

    /**
     * Id, mostly used for passing to R scripts
     */
    final int id

    RichnessEstimateType(int id) {
        this.id = id
    }
}
