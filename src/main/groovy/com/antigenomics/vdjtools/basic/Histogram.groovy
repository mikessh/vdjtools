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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.sample.Clonotype

abstract class Histogram {
    protected final boolean unweighted
    protected final int min, max, span

    protected final double[] innerHist
    protected final int[] lengths

    protected int count = 0
    protected double freq = 0

    public final String HEADER

    protected Histogram(boolean unweighted, int min, int max) {
        this.unweighted = unweighted

        this.min = min
        this.max = max

        this.span = max - min + 1
        this.innerHist = new double[span]
        this.lengths = (min..max) as double[]
        this.HEADER = lengths.collect().join("\t")
    }

    protected abstract int getValue(Clonotype clonotype)

    public int bin(Clonotype clonotype) {
        Math.min(max, Math.max(min, getValue(clonotype))) - min
    }

    /**
     * Update histogram with a single clonotype 
     * @param clonotype clonotype to add
     */
    public void add(Clonotype clonotype) {
        if (unweighted) {
            this.innerHist[bin(clonotype)]++
            freq++
        } else {
            this.innerHist[bin(clonotype)] += clonotype.freq
            freq += clonotype.freq
        }
        count++
    }

    /**
     * Update histogram with a set of clonotypes 
     * @param sample sample to add
     */
    public void addAll(Iterable<Clonotype> sample) {
        sample.each { Clonotype clonotype ->
            add(clonotype)
        }
    }

    /**
     * Combines two histograms
     * @param other histogram to update with
     */
    public void addAll(Histogram other) {
        if (other.unweighted != this.unweighted ||
                other.min != this.min || other.max != this.max)
            throw new Exception("Can't combine histograms of different type")

        this.count += other.count
        this.freq += other.freq
        other.innerHist.eachWithIndex { it, ind ->
            innerHist[ind] += it
        }
    }

    /**
     * Clears all spectratype counters
     */
    public void clear() {
        for (int i = 0; i < span; i++)
            this.innerHist[i] = 0
        count = 0
        freq = 0
    }

    /**
     * Normalize and get the histogram. Values correspond to specified clonotype feature, ranging from {@code min} to {@code max}
     * @return
     */
    public double[] getHistogram() {
        getHistogram(true)
    }

    /**
     * Normalize and get the histogram array. Values correspond to specified clonotype feature, ranging from {@code min} to {@code max}
     * @param normalized if true will normalize clonotype to total sum of 1.
     * @return
     */
    public double[] getHistogram(boolean normalized) {
        def histogram = new double[span]

        double _freq = (normalized && freq > 0) ? freq : 1.0

        for (int i = 0; i < span; i++)
            histogram[i] = this.innerHist[i] / _freq

        histogram
    }

    /**
     * Gets min CDR3 length of clonotype (lower boundary on spectratype length)
     * @return
     */
    int getMin() {
        min
    }

    /**
     * Gets max CDR3 length of clonotype (upper boundary on spectratype length)
     * @return
     */
    int getMax() {
        max
    }

    /**
     * Gets the number of bins in spectratype histogram, i.e. {@code max - min + 1}
     * @return
     */
    public int getSpan() {
        span
    }

    /**
     * Gets total frequency of summarized clonotypes
     * @return
     */
    public double getFreq() {
        freq
    }

    /**
     * Gets total read count of summarized clonotypes
     * @return
     */
    public int getCount() {
        count
    }

    /**
     * An array of CDR3 lengths that correspond to spectratype bins
     */
    public int[] getLengths() {
        return lengths
    }

    boolean isUnweighted() {
        unweighted
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        histogram.collect().join("\t")
    }
}
