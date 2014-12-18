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

import org.apache.commons.math3.util.CombinatoricsUtils
import sun.reflect.generics.reflectiveObjects.NotImplementedException

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
                        F0 +
                                F1 * (2 * F1 - 1) * (2 * F1 - 1) / 4 / (F2 + 1) / (F2 + 1) +
                                F1 * F1 * F2 * (F1 - 1) * (F1 - 1) / 4 / (F2 + 1) / (F2 + 1) / (F2 + 1) / (F2 + 1)
                ),
                frequencyTable.count,
                true, false, "chao1")
    }

    public Diversity chaoE(long extrapolateTo) {
        def mStar = extrapolateTo - n

        if (mStar < 0)
            throw new IllegalArgumentException("Should extrapolate farther than the total count of the sample")

        new Diversity(
                Sobs + F0 * (1.0 - Math.pow(1.0 - F1 / n / F0, mStar)),
                0, // todo: too lazy to differentiate right now
                extrapolateTo,
                true, false, "chao_e")
    }

    public Diversity chaoI(long interpolateTo) {
        if (n > Integer.MAX_VALUE)
            throw new NotImplementedException()

        if (interpolateTo < 0 || interpolateTo > n)
            throw new IllegalArgumentException("Should interpolate within the size of sample")

        double sum1 = 0, sum2 = 0, denom = CombinatoricsUtils.binomialCoefficientLog((int) n, (int) interpolateTo)

        frequencyTable.bins.each { FrequencyTable.FrequencyTableBin bin ->
            def k = bin.count,
                f = bin.diversity

            if (k <= n - interpolateTo) {
                def alpha = Math.exp(CombinatoricsUtils.binomialCoefficientLog((int) n - k, (int) interpolateTo) - denom)
                sum1 += f * alpha
                sum2 += f * (1 - alpha) * (1 - alpha)
            }
        }

        double Sind = Sobs - sum1

        new Diversity(
                Sind,
                Math.sqrt(sum2 - Sind * Sind / (Sobs + F0)),
                interpolateTo,
                false, false, "chao_i")
    }
}
