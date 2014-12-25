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
 * Last modified on 15.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.Sample

class Rarefaction {
    private final FrequencyTable frequencyTable
    private final ChaoEstimator chaoEstimator
    private final long n


    public Diversity getAt(long coord) {
        coord > n ? chaoEstimator.chaoE(coord) : chaoEstimator.chaoI(coord)
    }

    public Rarefaction(Sample sample, IntersectionType intersectionType) {
        this.frequencyTable = new FrequencyTable(sample, intersectionType)
        this.n = frequencyTable.count
        this.chaoEstimator = new ChaoEstimator(frequencyTable)
    }

    public Rarefaction(Sample sample) {
        this(sample, IntersectionType.Strict)
    }

    public ArrayList<RarefactionPoint> interpolate() {
        build(n)
    }

    public ArrayList<RarefactionPoint> extrapolate(long to) {
        build(n + 1, to)
    }

    public ArrayList<RarefactionPoint> interpolate(int numberOfPoints) {
        build(n, numberOfPoints)
    }

    public ArrayList<RarefactionPoint> extrapolate(long to, int numberOfPoints) {
        build(n + 1, to, numberOfPoints)
    }

    public ArrayList<RarefactionPoint> build(long to) {
        build(0L, to)
    }

    public ArrayList<RarefactionPoint> build(long from, long to) {
        build(from, to, Math.min(101, (int) n))
    }

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

        // add the last point (more robust)
        rarefactionCurve.add(new RarefactionPoint(this[to]))

        rarefactionCurve.sort { it.x }
    }

    class RarefactionPoint {
        public final double x, mean, ciU, ciL
        public final DiversityType diversityType

        RarefactionPoint(Diversity diversity) {
            this.diversityType = diversity.type
            this.x = diversityType == DiversityType.TotalDiversityLowerBoundEstimate ? 1e20 : diversity.n
            this.mean = diversity.mean
            this.ciL = mean - 1.96 * diversity.std
            this.ciU = mean + diversity.std
        }

        public static final String HEADER = "x\tmean\tciL\tciU\ttype"

        @Override
        public String toString() {
            [x, mean, ciL, ciU, diversityType.id].join("\t")
        }
    }
}
