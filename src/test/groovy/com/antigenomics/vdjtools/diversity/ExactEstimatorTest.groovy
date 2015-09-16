/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
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

                case "shannonWienerIndex":
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

        def trials = 300

        //println DiversityEstimator.HEADER
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
            assert rate >= 0.7
        }
    }
}
