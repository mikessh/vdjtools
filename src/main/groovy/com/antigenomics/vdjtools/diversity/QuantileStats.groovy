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
 * Last modified on 3.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.ClonotypeContainer
import com.google.common.util.concurrent.AtomicDouble
import com.google.common.util.concurrent.AtomicDoubleArray


class QuantileStats {
    private final int numberOfQuantiles
    private final AtomicDoubleArray quantileFreqs
    private final AtomicDouble highOrderFreq = new AtomicDouble(),
                               doubletonFreq = new AtomicDouble(),
                               singletonFreq = new AtomicDouble()
    private final double totalFreq

    public QuantileStats(ClonotypeContainer clonotypeContainer, int numberOfQuantiles) {
        this.numberOfQuantiles = numberOfQuantiles
        this.quantileFreqs = new AtomicDoubleArray(numberOfQuantiles)

        if (!clonotypeContainer.isSorted())
            throw new Exception("Clonotype container should be sorted to be used as input for this statistic")

        this.totalFreq = clonotypeContainer.freq

        update(clonotypeContainer)
    }

    public QuantileStats(ClonotypeContainer clonotypeContainer) {
        this(clonotypeContainer, 5)
    }

    private void update(ClonotypeContainer clonotypeContainer) {
        int n = clonotypeContainer.diversity, m = n

        for (int i = n - 1; i >= 0; i--) {
            def clonotype = clonotypeContainer[i]
            def count = clonotype.getCount(),
                freq = clonotype.getFreq()
            boolean highOrderFlag = false
            switch (count) {
                case 1:
                    singletonFreq.addAndGet(freq)
                    break
                case 2:
                    doubletonFreq.addAndGet(freq)
                    break
                default:
                    highOrderFlag = true
                    break
            }
            if (highOrderFlag) {
                m = i
                break
            }
        }

        for (int i = 0; i <= m; i++) {
            int q = (i * numberOfQuantiles) / (m + 1)
            def clonotype = clonotypeContainer[i]
            def freq = clonotype.getFreq()
            highOrderFreq.addAndGet(freq)
            quantileFreqs.addAndGet(q, freq)
        }
    }

    public int getNumberOfQuantiles() {
        return numberOfQuantiles
    }

    public double getQuantileFrequency(int quantile) {
        if (quantile < 0 || quantile >= numberOfQuantiles)
            throw new IndexOutOfBoundsException()
        quantileFreqs.get(quantile) / totalFreq
    }

    public double getSingletonFreq() {
        singletonFreq.get() / totalFreq
    }

    public double getDoubletonFreq() {
        doubletonFreq.get() / totalFreq
    }

    public double getHighOrderFreq() {
        highOrderFreq.get() / totalFreq
    }

    public final String HEADER = "singleton\tdoubleton\t" +
            (1..numberOfQuantiles).collect { "high_order_Q$it" }.join("\t")

    @Override
    String toString() {
        [singletonFreq, doubletonFreq, (0..<numberOfQuantiles).collect { getQuantileFrequency(it) }].flatten().join("\t")
    }
}
