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
 * Last modified on 31.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.intersection

/**
 * Normalization type that is recommended for a given IntersectMetric
 * Log: log10(x+1e-9)
 * Correlation: (1-x)/2
 * None: x 
 */
public enum IntersectMetricNormalization {
    Log(0), Correlation(1), None(2)

    /**
     * Normalization type ID, used for passing to R scripts as argument
     */
    public final int id

    private IntersectMetricNormalization(int id) {
        this.id = id
    }

    public double normalize(double x) {
        switch (this) {
            case Log:
                return Math.log10(x + 1.0)
            case Correlation:
                return (1.0 - x) / 2
        }
        x
    }
}