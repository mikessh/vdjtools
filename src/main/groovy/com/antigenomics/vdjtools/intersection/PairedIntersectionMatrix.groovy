/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.sample.SampleCollection

class PairedIntersectionMatrix {
    public final IntersectionUtil intersectionUtil
    private final List<PairedIntersection> pairedIntersections
    private final PairedIntersection[][] intersectionMatrix
    private final SampleCollection parentCollection
    private final int n

    PairedIntersectionMatrix(SampleCollection parentCollection,
                             List<PairedIntersection> pairedIntersections,
                             IntersectionUtil intersectionUtil) {
        this.intersectionUtil = intersectionUtil
        this.pairedIntersections = pairedIntersections
        this.parentCollection = parentCollection
        this.n = parentCollection.size()
        this.intersectionMatrix = new PairedIntersection[n][n]
        pairedIntersections.each {
            intersectionMatrix[it.parent.i][it.parent.j] = it
        }
    }

    PairedIntersection getAt(int i, int j) {
        if (i >= j || j >= n)
            throw new IllegalArgumentException("i < j < $n required")
        intersectionMatrix[i][j]
    }

    /**
     * Builds a symmetric intersection matrix using provided metric.
     * Diagonal elements are masked with NaNs.
     * @param metric metric to characterize the degree of intersection between sample pair
     * @return symmetric intersection matrix
     */
    double[][] buildDistanceMatrix(IntersectMetric metric) {
        def matrix = new double[n][n]
        for (int i = 0; i < n; i++) {
            matrix[i][i] = Double.NaN

            for (int j = i + 1; j < n; j++) {
                matrix[i][j] = metric.value(intersectionMatrix[i][j])
                matrix[j][i] = matrix[i][j]
            }
        }

        matrix
    }

    /**
     * Builds an frequency matrix, where (i,j) and (j,i) elements corresponds to
     * frequency of (i,j) overlapping clonotypes in i-th and j-th samples respectively.
     * The resulting matrix is non-symmetric, and thus formally can't be used as metric,
     * see buildIntersectMatrix.
     * Diagonal elements are masked with NaNs.
     * @return non-symmetric intersection matrix
     */
    double[][] buildFrequencyMatrix() {
        def matrix = new double[n][n]

        for (int i = 0; i < n; i++) {
            matrix[i][i] = Double.NaN

            for (int j = i + 1; j < n; j++) {
                matrix[i][j] = intersectionMatrix[i][j].freq12 // intersectionMatrix[i][j].freq12e
                matrix[j][i] = intersectionMatrix[i][j].freq21 // intersectionMatrix[i][j].freq21e
            }
        }

        matrix
    }

    int size() {
        n
    }

    @Override
    String toString() {
        "Paired intersection of $n samples"
    }

    void print(PrintWriter pw, boolean header = true) {
        if (header)
            pw.println("#" +
                    [PairedIntersection.HEADER,
                     parentCollection.metadataTable.columnHeader1,
                     parentCollection.metadataTable.columnHeader2].flatten().join("\t").trim())

        pairedIntersections.each { PairedIntersection pairedIntersection ->
            pw.println([pairedIntersection.toString(),
                        pairedIntersection.sample1.sampleMetadata.toString(),
                        pairedIntersection.sample2.sampleMetadata.toString()].join("\t").trim())
        }
    }
}
