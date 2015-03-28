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

package com.antigenomics.vdjtools.overlap

/**
 * Normalization type that is recommended for a given IntersectMetric,
 * should transform overlap metric value to {@code [0, +inf)} scale.
 */
public enum OverlapMetricNormalization {
    /**
     * Negative logarithm normalization {@code -log10(x + 1e-9)}.
     */
    /*   */ NegLog(0),
    /**
     * Correlation metrics, normalized as {@code ( 1 - x ) /2}.
     */
            R(1),
    /**
     * Metrics for which normalization is not required.
     */
            None(2)

    /**
     * Normalization type ID, used for passing to R scripts as argument.
     */
    public final int id

    private OverlapMetricNormalization(int id) {
        this.id = id
    }

    /**
     * Normalize raw overlap metric accordingly.
     * @param x overlap metric value.
     * @return normalized value.
     */
    public double normalize(double x) {
        switch (this) {
            case NegLog:
                return -Math.log10(x + 1e-9)
            case R:
                return (1.0 - x) / 2
        }
        x
    }
}