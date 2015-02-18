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
 * Last modified on 16.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.ClonotypeContainer
import com.google.common.util.concurrent.AtomicDouble
import com.google.common.util.concurrent.AtomicDoubleArray


/**
 * A class that computes summary statistics of repertoire clonality divided into several levels:
 * 1. singletons (encountered once), doubletons (encountered twice) and high-order clonotypes - 
 *    those are base quantities to estimate the lower bound on total repertoire diversity
 * 2. cumulative frequency for several quantiles (e.g. top 25%, next 25%, ...) of high-order clonotypes
 * 3. details for top N clonotypes 
 */
public class QuantileStats {
    private final int numberOfQuantiles
    private final AtomicDoubleArray quantileFreqs
    private final AtomicDouble highOrderFreq = new AtomicDouble(),
                               doubletonFreq = new AtomicDouble(),
                               singletonFreq = new AtomicDouble()
    private final double totalFreq

    /**
     * Summarizes quantile statisitcs for a given sample 
     * @param clonotypeContainer a set of clonotypes
     * @param numberOfQuantiles number of quantiles for 2nd level of detalizaiton
     */
    public QuantileStats(ClonotypeContainer clonotypeContainer, int numberOfQuantiles) {
        this.numberOfQuantiles = numberOfQuantiles
        this.quantileFreqs = new AtomicDoubleArray(numberOfQuantiles)

        if (!clonotypeContainer.isSorted())
            throw new Exception("Clonotype container should be sorted to be used as input for this statistic")

        this.totalFreq = clonotypeContainer.count / (double) clonotypeContainer.diversity

        update(clonotypeContainer)
    }

    /**
     * * Summarizes quantile statisitcs for a given sample 
     * @param clonotypeContainer a set of clonotypes
     */
    public QuantileStats(ClonotypeContainer clonotypeContainer) {
        this(clonotypeContainer, 5)
    }

    /**
     * Internal - adds more clonotyps to stats
     */
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

    /**
     * Gets the number of 2nd level summary quantiles 
     * @return
     */
    public int getNumberOfQuantiles() {
        return numberOfQuantiles
    }

    /**
     * Gets frequency for a given quantile
     * @param quantile quantile index, should be less than {@code numberOfQuantiles} and 
     *                 greater or equal than {@code 0}
     * @return
     * @throws IndexOutOfBoundsException
     */
    public double getQuantileFrequency(int quantile) {
        if (quantile < 0 || quantile >= numberOfQuantiles)
            throw new IndexOutOfBoundsException()
        quantileFreqs.get(quantile) / totalFreq
    }

    /**
     * Gets the frequency of singletons, i.e. clonotypes represented by a single read 
     * @return
     */
    public double getSingletonFreq() {
        singletonFreq.get() / totalFreq
    }

    /**
     * Gets the frequency of doubletons, i.e. clonotypes represented by two reads 
     * @return
     */
    public double getDoubletonFreq() {
        doubletonFreq.get() / totalFreq
    }

    /**
     * Gets the frequency of high order clonotypes, i.e. clonotypes represented by more than two reads
     * @return
     */
    public double getHighOrderFreq() {
        highOrderFreq.get() / totalFreq
    }

    /**
     * Header string, used for tabular output
     */
    public static final String HEADER = "type\tname\tvalue"

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        ["set\t3+\t$highOrderFreq",
         "set\t2\t$doubletonFreq",
         "set\t1\t$singletonFreq",
         (0..<numberOfQuantiles).collect { "quantile\tQ${it + 1}\t${getQuantileFrequency(it)}" }].flatten().join("\n")
    }
}
