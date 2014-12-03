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

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.ClonotypeContainer
import com.google.common.util.concurrent.AtomicDouble
import com.google.common.util.concurrent.AtomicDoubleArray


class QuantileStats {
    private final int numberOfQuantiles
    private final AtomicDoubleArray quantileFreqs
    private final AtomicDouble sumFreq = new AtomicDouble(0)

    public QuantileStats(ClonotypeContainer clonotypeContainer, int numberOfQuantiles) {
        this.numberOfQuantiles = numberOfQuantiles
        this.quantileFreqs = new AtomicDoubleArray(numberOfQuantiles)

        if (!clonotypeContainer.isSorted())
            throw new Exception("Clonotype container should be sorted to be used as input for this statistic")

        update(clonotypeContainer)
    }

    public QuantileStats(ClonotypeContainer clonotypeContainer) {
        this(clonotypeContainer, 5)
    }

    private void update(ClonotypeContainer clonotypeContainer) {
        final int n = clonotypeContainer.diversity

        clonotypeContainer.eachWithIndex { Clonotype clonotype, int ind ->
            int q = (ind * numberOfQuantiles) / n
            double freq = clonotype.freq
            quantileFreqs.addAndGet(q, freq)
            sumFreq.addAndGet(freq)
        }
    }

    public int getNumberOfQuantiles() {
        return numberOfQuantiles
    }

    public double getFrequency(int quantile) {
        if (quantile < 0 || quantile >= numberOfQuantiles)
            throw new IndexOutOfBoundsException()
        quantileFreqs.get(quantile) / sumFreq.get()
    }

    public final String HEADER = (1..numberOfQuantiles).collect { "Q$it" }.join("\t")

    @Override
    String toString() {
        (0..<numberOfQuantiles).collect { getFrequency(it) }.join("\t")
    }
}
