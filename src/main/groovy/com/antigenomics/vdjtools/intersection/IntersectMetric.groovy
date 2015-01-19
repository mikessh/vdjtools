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

package com.antigenomics.vdjtools.intersection

import static com.antigenomics.vdjtools.intersection.IntersectMetricNormalization.*

/**
 * An enum that defines intersection metric, a function that characterizes the extent of overlap between a pair of samples.
 */
public enum IntersectMetric {
    /**
     *
     */
    Correlation("R", Correlation),
    /**
     *
     */
            Diversity("D", NegLog),
    Frequency("F", NegLog),
    Frequency2("F2", NegLog),
    vJSD("vJSD", None),
    vjJSD("vjJSD", None),
    vj2JSD("vj2JSD", None),
    sJSD("sJSD", None)

    public final String shortName
    public final IntersectMetricNormalization normalization

    /**
     * Defines a new intersection metric 
     * @param shortName short name
     * @param normalization normalization type
     */
    public IntersectMetric(String shortName, IntersectMetricNormalization normalization) {
        this.shortName = shortName
        this.normalization = normalization
    }

    /**
     * Gets {@code IntersectMetric} by short name
     * @param shortName short name
     * @return
     */
    public static IntersectMetric getByShortName(String name) {
        name = name.toUpperCase()
        values().find { it.shortName.toUpperCase() == name }
    }

    /**
     * A list of existing {@code IntersectMetric} short names
     */
    public static String allowedNames = values().collect { it.shortName }.join(",")
}