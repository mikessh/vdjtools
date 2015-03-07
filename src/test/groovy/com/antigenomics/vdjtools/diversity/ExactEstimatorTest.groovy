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
 * Last modified on 6.3.2015 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import org.junit.Test

class ExactEstimatorTest {

    private static boolean check(String estimateName,
                                 DiversityEstimate estimate,
                                 FrequencyTableGenerator tableGenerator) {
        def lower = estimate.mean - estimate.std,
            upper = estimate.mean + estimate.std

        if (estimate instanceof DiversityIndex) {
            switch (estimateName) {
                case "d50Index":
                    assert lower < 1.001
                    return upper > 0.9

                case "shannonWeinerIndex":
                    assert lower < 1.001 * tableGenerator.observedSpecies
                    return upper > 0.7 * tableGenerator.observedSpecies

                case "inverseSimpsonIndex":
                    assert lower < 1.001 * tableGenerator.observedSpecies
                    return upper > 0.4 * tableGenerator.observedSpecies
            }
        } else if (estimate instanceof SpeciesRichness) {
            switch (estimateName) {
                case "observedDiversity":
                    assert upper > 0.999 * tableGenerator.observedSpecies &&
                            lower < 1.001 * tableGenerator.observedSpecies
                    return true

                case "chaoE":
                    return upper > 0.5 * tableGenerator.numberOfSpecies &&
                            lower < 1.001 * tableGenerator.numberOfSpecies

                case "efronThisted":
                    return upper > 0.5 * tableGenerator.numberOfSpecies &&
                            lower < 1.001 * tableGenerator.numberOfSpecies

                case "chao1":
                    return upper > 0.5 * tableGenerator.numberOfSpecies &&
                            lower < 1.001 * tableGenerator.numberOfSpecies
            }
        }
    }

    @Test
    public void randomTest() {
        def tableGenerator = new FrequencyTableGenerator()

        def goodCounter = new HashMap<String, Integer>()
        DiversityEstimator.ESTIMATE_NAMES.each { goodCounter.put(it, 0) }

        def trials = 100

        println DiversityEstimator.HEADER
        for (int i = 0; i < trials; i++) {
            def table = tableGenerator.create()

            def diversityEstimates = new ExactEstimator(table, table.count * 10)

            diversityEstimates.computeAll().each {
                def name = it.key, estimate = it.value, good = check(name, estimate, tableGenerator)

                //println it.toString() + "\t" + good

                if (good)
                    goodCounter.put(name, (goodCounter[name] ?: 0) + 1)

            }

            //println diversityEstimates
        }

        goodCounter.each {
            def rate = it.value / (double) trials
            println "Good estimate rate for $it.key is $rate"
            assert rate >= 0.8
        }
    }
}
