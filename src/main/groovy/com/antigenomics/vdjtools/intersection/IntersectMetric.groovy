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



package com.antigenomics.vdjtools.intersection

import static com.antigenomics.vdjtools.intersection.IntersectMetricNormalization.*

/**
 * An enum that defines intersection metric, a function that characterizes the extent of overlap between a pair of samples.
 */
public enum IntersectMetric {
    /**
     * Correlation between sample frequencies of _overlapping_ clonotypes
     */
    /*   */ Correlation("R", R),
    /**
     * Ratio of observed to expected numbers of unique overlapping clonotypes, {@code div12 / div1 / div2}
     */
            Diversity("D", NegLog),
    /**
     * Geometric mean of sums of frequencies of overlapping clonotypes, {@code sqrt ( freq12 * freq21 )}
     */
            Frequency("F", NegLog),
    /**
     * Sum of geometric means of frequencies of overlapping clonotypes, {@code sum ( sqrt ( freq12(i) * freq21(i) ) ), i = 1..div12}
     */
            Frequency2("F2", NegLog),
    /**
     * Jensen-Shannon divergence between Variable segment usage vectors
     */
            vJSD("vJSD", None),
    /**
     * Jensen-Shannon divergence between concatenated Variable and Joining segment usage vectors
     */
            vjJSD("vjJSD", None),
    /**
     * Jensen-Shannon divergence between flattened Variable-Joining segment pairing matrices
     */
            vj2JSD("vj2JSD", None),
    /**
     * Jensen-Shannon divergence between spectratypes
     */
            sJSD("sJSD", None),
    /**
     * Jaccard index
     */
            Jaccard("Jaccard", None),
    /**
     * MorisitaHorn index
     */
            MorisitaHorn("MorisitaHorn", None)

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