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
 * Last modified on 7.3.2015 by mikesh
 */

package com.antigenomics.vdjtools.diversity

/**
 * Specifies how the diversity estimate was computed
 */
enum EstimationMethod {
    /**
     * Exact method was used to compute mean and standard deviation 
     */
            Exact("exact"),
    /**
     * Re-sampling was used to compute mean and standard deviation 
     */
            Resampled("resampled")
    
    private final String name

    EstimationMethod(String name) {
        this.name = name
    }

    String getName() {
         name
    }

    @Override
    String toString() {
        name
    }
}