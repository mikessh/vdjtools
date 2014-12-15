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


class RarefactionCurve {
    private final int[][] y
    private final int[] mean, std, x
    private final int numberOfResamples, numberOfPoints
    public final static double WINDOW_SIZE_FACTOR = 0.05

    RarefactionCurve(int[] x, int[][] y) {
        this.y = y
        this.x = x
        this.numberOfResamples = y[0].length
        this.numberOfPoints = x.length

        this.mean = new int[numberOfPoints]
        this.std = new int[numberOfPoints]

        def mean = new int[numberOfPoints],
            std = new int[numberOfPoints]

        int windowHalfSz = Math.max(1, (int)(numberOfPoints * WINDOW_SIZE_FACTOR)), windowSz = 2 * windowHalfSz + 1

        for (int i = 0; i < numberOfPoints; i++) {
            double yy = 0, yy2 = 0
            for (int j = 0; j < numberOfResamples; j++) {
                yy += y[i][j]
                yy2 += y[i][j] * y[i][j]
            }
            mean[i] = yy / numberOfResamples
            std[i] = Math.sqrt((yy2 - yy * yy / numberOfResamples) / (numberOfResamples - 1))
        }

        // sliding mean smoothing for now / GPR for later
        for (int i = 1; i < numberOfPoints - 1; i++) {
            int mSum = 0, sSum = 0
            for (int j = i - windowHalfSz; j <= i + windowHalfSz; j++) {
                if (j < 0) {
                    mSum += 0
                    sSum += 0
                } else if (j >= numberOfPoints) {
                    mSum += mean[numberOfPoints - 1]
                    sSum += 0
                } else {
                    mSum += mean[j]
                    sSum += std[j]
                }
            }
            this.mean[i] = mSum / windowSz
            this.std[i] = sSum / windowSz
        }
        this.mean[numberOfPoints - 1] = mean[numberOfPoints - 1]
    }

    int getNumberOfPoints() {
        numberOfPoints
    }

    RarefactionPoint getAt(int i) {
        int m = mean[i], s = std[i]
        new RarefactionPoint(x[i], m, m - s, m + s, i == numberOfPoints - 1)
    }

    class RarefactionPoint {
        public final int x, mean, ciU, ciL
        public final boolean last

        RarefactionPoint(int x, int mean, int ciL, int ciU, boolean last) {
            this.x = x
            this.mean = mean
            this.ciL = ciL
            this.ciU = ciU
            this.last = last
        }

        public static final String HEADER = "x\tmean\tciL\tciU\tlast"

        @Override
        public String toString() {
            [x, mean, ciL, ciU, last ? 1 : 0].join("\t")
        }
    }
}
