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
 * Last modified on 18.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.util.ExecUtil
import com.google.common.util.concurrent.AtomicDouble
import groovyx.gpars.GParsPool
import org.apache.commons.math3.util.CombinatoricsUtils
import sun.reflect.generics.reflectiveObjects.NotImplementedException

import static com.antigenomics.vdjtools.diversity.DiversityType.*

class ChaoEstimator {
    private final FrequencyTable frequencyTable
    private final long n
    private final double F0, Sobs, F1, F2

    public ChaoEstimator(FrequencyTable frequencyTable) {
        this.frequencyTable = frequencyTable
        this.n = frequencyTable.count
        this.Sobs = frequencyTable.diversity
        this.F1 = frequencyTable.singletons
        this.F2 = frequencyTable.doubletons
        this.F0 = F1 * (F1 - 1) / 2 / (F2 + 1)
    }

    public Diversity chao1() {
        new Diversity(
                Sobs + F0,
                Math.sqrt(
                        F0 + F1 * (2 * F1 - 1) * (2 * F1 - 1) / 4 / (F2 + 1) / (F2 + 1) +
                                F1 * F1 * F2 * (F1 - 1) * (F1 - 1) / 4 / (F2 + 1) / (F2 + 1) / (F2 + 1) / (F2 + 1)
                ),
                frequencyTable.count,
                TotalDiversityLowerBoundEstimate, false, "chao1")
    }

    public Diversity chaoE(long extrapolateTo) {
        double mStar = extrapolateTo - n

        if (mStar < 0)
            throw new IllegalArgumentException("Should extrapolate farther than the total count of the sample")

        // here is the derivative to compute variance
        double dF0dF1 = (2 * F1 - 1) / 2 / (F2 + 1),
               dF0dF2 = -F1 * (F1 - 1) / 2 / (F2 + 1) / (F2 + 1)
        double brackets = 1.0 - Math.pow(1.0 - F1 / n / F0, mStar),
               dBrackets = -mStar * Math.pow(1.0 - F1 / n / F0, mStar - 1),
               dBracketsdF1 = dBrackets * (-1 / n / F0 + F1 / n / F0 / F0 * dF0dF1),
               dBracketsdF2 = dBrackets * (F1 / n / F0 / F0 * dF0dF2)
        double dSdF1 = dF0dF1 * brackets + F0 * dBracketsdF1,
               dSdF2 = dF0dF2 * brackets + F0 * dBracketsdF2,
               cov11 = F1 * (1 - F1 / (Sobs + F0)), cov22 = F2 * (1 - F2 / (Sobs + F0)),
               cov12 = -F1 * F2 / (Sobs + F0)

        new Diversity(
                Sobs + F0 * brackets,
                Math.sqrt(dSdF1 * dSdF1 * cov11 + dSdF2 * dSdF2 * cov22 + 2 * dSdF1 * dSdF2 * cov12),
                extrapolateTo,
                extrapolateTo == n ? Observed : Extrapolated, false, "chao_e")
    }

    public Diversity chaoI(long interpolateTo) {
        if (n > Integer.MAX_VALUE)
            throw new NotImplementedException()

        if (interpolateTo < 0 || interpolateTo > n)
            throw new IllegalArgumentException("Should interpolate within the size of sample")


        double denom = CombinatoricsUtils.binomialCoefficientLog((int) n, (int) interpolateTo)
        def sum1 = new AtomicDouble(), sum2 = new AtomicDouble()

        GParsPool.withPool ExecUtil.THREADS, { // this is quite time consuming
            frequencyTable.bins.eachParallel { FrequencyTable.FrequencyTableBin bin ->
                def k = bin.count,
                    f = bin.diversity

                if (k <= n - interpolateTo) {
                    def alpha = Math.exp(CombinatoricsUtils.binomialCoefficientLog((int) n - k,
                            (int) interpolateTo) - denom)
                    sum1.addAndGet(f * alpha)
                    sum2.addAndGet(f * (1 - alpha) * (1 - alpha))
                }
            }
        }

        sum1 = sum1.get()
        sum2 = sum2.get()

        double Sind = Sobs - sum1

        new Diversity(
                Sind,
                Math.sqrt(sum2 - Sind * Sind / (Sobs + F0)),
                interpolateTo,
                interpolateTo == n ? Observed : Interpolated, false, "chao_i")
    }
}
