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

import org.apache.commons.math3.distribution.ZipfDistribution

class FrequencyTableGenerator {
    int numberOfSpecies = 500, observedSpecies = 100
    double exponent = 3.0

    FrequencyTableGenerator() {


    }

    FrequencyTableGenerator(int numberOfSpecies, double exponent) {
        this.numberOfSpecies = numberOfSpecies
        this.exponent = exponent
    }

    public FrequencyTable create() {
        def distr = new ZipfDistribution(numberOfSpecies, exponent)

        def cache = new HashMap<Long, Long>()

        observedSpecies.times {
            def count = distr.sample()
            cache.put(count, (cache[count] ?: 0) + 1)
        }

        new FrequencyTable(cache)
    }
}
