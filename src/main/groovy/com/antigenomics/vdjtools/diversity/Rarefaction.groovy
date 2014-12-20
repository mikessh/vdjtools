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

    public ArrayList<RarefactionPoint> build(long extrapolateTo) {
        build(extrapolateTo, Math.min(100, (int) n))
    }

    public ArrayList<RarefactionPoint> build(long extrapolateTo, int numberOfPoints) {
        def rarefactionCurve = new ArrayList<RarefactionPoint>()

        int step = extrapolateTo / (numberOfPoints - 1)
        boolean hasExact = false

        for (int i = 0; i < numberOfPoints; i++) {
            long m = i * step

            if (m == n)
                hasExact = true

            rarefactionCurve.add(new RarefactionPoint(this[m]))
        }

        if (!hasExact)
            rarefactionCurve.add(new RarefactionPoint(chaoEstimator.chaoI(n)))

        rarefactionCurve.sort { it.x }
    }

    class RarefactionPoint {
        public final double x, mean, ciU, ciL
        public final DiversityType diversityType

        RarefactionPoint(Diversity diversity) {
            this.diversityType = diversity.type
            this.x = diversityType == DiversityType.TotalDiversityLowerBoundEstimate ? 1e20 : diversity.n
            this.mean = diversity.mean
            this.ciL = mean - diversity.std
            this.ciU = mean + diversity.std
        }

        public static final String HEADER = "x\tmean\tciL\tciU\ttype"

        @Override
        public String toString() {
            [x, mean, ciL, ciU, diversityType.id].join("\t")
        }
    }
}
