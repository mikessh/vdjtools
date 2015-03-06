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
 * Last modified on 19.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.Sample

/**
 * A class that is used to compute a rarefaction curve based on {@code ChaoEstimator}. 
 * Rarefaction curve shows the relation between sample diversity and sample size.
 * Implementation is based on Colwell et al 2012 paper.
 * 
 * @url http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/References/ColwellEtAl2012.pdf
 */
public class Rarefaction {
    private final FrequencyTable frequencyTable
    private final ChaoEstimator chaoEstimator
    private final long n

    /**
     * Computes the interpolated/exact/extrapolated diversity estimate for a given sample size
     * @param coord sample size
     * @return
     */
    public Diversity getAt(long coord) {
        coord > n ? chaoEstimator.chaoE(coord) : chaoEstimator.chaoI(coord)
    }

    /**
     * Creates an instance that will perform rarefaction for a given sample  
     * @param sample sample to analyze
     * @param intersectionType {@code IntersectionType} that will be used to collapse clonotypes when computing a {@code FrequencyTable}
     */
    public Rarefaction(Sample sample, IntersectionType intersectionType) {
        this.frequencyTable = new FrequencyTable(sample, intersectionType)
        this.n = frequencyTable.count
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Creates an instance that will perform rarefaction for a given frequency table 
     * @param frequencyTable {@code FrequencyTable} summarizing the number of clonotypes met once, twice, ... in a given clonotype set
     */
    public Rarefaction(FrequencyTable frequencyTable) {
        this.frequencyTable = frequencyTable
        this.n = frequencyTable.count
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    /**
     * Creates an instance that will perform rarefaction for a given sample. 
     * Will use {@code IntersectionType.strict} to build a frequency table 
     * @param sample sample to analyze
     */
    public Rarefaction(Sample sample) {
        this(sample, IntersectionType.Strict)
    }

    /**
     * Builds a rarefaction curve interpolation up to given sample size 
     * @return
     */
    public ArrayList<RarefactionPoint> interpolate() {
        build(n)
    }

    /**
     * Extrapolates rarefaction curve from the size of a given sample up to a specified size
     * @param to where to extrapolate 
     * @return
     * @throws IllegalArgumentException if {@code to} is less than sample size + 1 
     */
    public ArrayList<RarefactionPoint> extrapolate(long to) {
        build(n + 1, to)
    }

    /**
     * Builds a rarefaction curve interpolation up to given sample size
     * @param numberOfPoints number of size steps
     * @return
     */
    public ArrayList<RarefactionPoint> interpolate(int numberOfPoints) {
        build(n, numberOfPoints)
    }

    /**
     * Extrapolates rarefaction curve starting from the size of a given sample up to a specified size
     * @param to where to extrapolate 
     * @param numberOfPoints number of points in rarefaction curve
     * @return
     * @throws IllegalArgumentException if {@code to} is less than sample size + 1 
     */
    public ArrayList<RarefactionPoint> extrapolate(long to, int numberOfPoints) {
        build(n + 1, to, numberOfPoints)
    }

    /**
     * Builds a rarefaction curve from sample size {@code 0} up to sample size {@code to}
     * @param to sample size for last rarefaction point
     * @return
     */
    public ArrayList<RarefactionPoint> build(long to) {
        build(0L, to)
    }

    /**
     * Builds a rarefaction curve starting from sample size {@code from} up to sample size {@code to}
     * @param from sample size for first rarefaction point
     * @param to sample size for last rarefaction point
     * @return
     */
    public ArrayList<RarefactionPoint> build(long from, long to) {
        build(from, to, Math.min(101, (int) n))
    }

    /**
     * Builds a rarefaction curve starting from sample size {@code from} up to sample size {@code to}
     * @param from sample size for first rarefaction point
     * @param to sample size for last rarefaction point
     * @param numberOfPoints number of points in rarefaction curve
     * @return
     */
    public ArrayList<RarefactionPoint> build(long from, long to, int numberOfPoints) {
        if (from > to)
            throw new IllegalArgumentException("from should be less than to")
        def rarefactionCurve = new ArrayList<RarefactionPoint>()

        double step = (to - from) / (double) (numberOfPoints - 1)
        boolean hasExact = false

        for (int i = 0; i < numberOfPoints - 1; i++) {
            long m = from + i * step

            if (m == n)
                hasExact = true

            rarefactionCurve.add(new RarefactionPoint(this[m]))
        }

        if (!hasExact)
            rarefactionCurve.add(new RarefactionPoint(chaoEstimator.chaoI(n)))

        // add the last point (more robust & label placing in plotting)
        rarefactionCurve.add(new RarefactionPoint(chaoEstimator.chaoE(to)))

        rarefactionCurve.sort { it.x }
    }

    /**
     * Holds summary for rarefied diversity estimate 
     */
    public class RarefactionPoint {
        private final double x, mean, ciU, ciL
        private final DiversityType diversityType

        public RarefactionPoint(Diversity diversity) {
            this.diversityType = diversity.type
            this.x = diversityType == DiversityType.TotalDiversityLowerBoundEstimate ? 1e20 : diversity.n
            this.mean = diversity.mean
            this.ciL = mean - 1.96 * diversity.std
            this.ciU = mean + 1.96 * diversity.std
        }

        /**
         * Gets the coordinate, i.e. sample size 
         * @return
         */
        public double getX() {
            x
        }

        /**
         * Gets the mean of rarefied diversity estimate 
         * @return
         */
        public double getMean() {
            mean
        }

        /**
         * Gets the upper bound of 95% confidence interval of rarefied diversity estimate
         * @return
         */
        public double getCiU() {
            ciU
        }

        /**
         * Gets the lower bound of 95% confidence interval of rarefied diversity estimate
         * @return
         */
        public double getCiL() {
            ciL
        }

        /**
         * Gets the diversity estimate type: interpolated, exact or extrapolated
         * @return
         */
        public DiversityType getDiversityType() {
            diversityType
        }

        /**
         * Header string, used for tabular output
         */
        public static final String HEADER = "x\tmean\tciL\tciU\ttype"

        /**
         * Plain text row for tabular output
         */
        @Override
        public String toString() {
            [x, mean, ciL, ciU, diversityType.id].join("\t")
        }
    }
}
