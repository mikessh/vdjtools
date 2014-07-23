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

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.timecourse.DynamicClonotype
import com.antigenomics.vdjtools.timecourse.TimeCourse

// todo: summary

class SequentialIntersection {
    private final Sample[] samples
    private final PairedIntersection[] pairedIntersections
    private final PairedIntersection[][] pairedIntersectionsMat
    private final IntersectionUtil intersectionUtil
    private boolean allIntersectionsBuilt = false

    SequentialIntersection(Sample[] samples, IntersectionType intersectionType) {
        if (samples.size() < 3)
            throw new IllegalArgumentException("More than 2 samples should be provided")

        this.samples = samples
        this.intersectionUtil = new IntersectionUtil(intersectionType)

        this.pairedIntersections = (0..(samples.length - 2)).collect { int i ->
            intersectionUtil.generatePairedIntersection(samples[i], samples[i + 1])
        } as PairedIntersection[]

        this.pairedIntersectionsMat = new PairedIntersection[samples.size()][samples.size()]
    }

    double[][] buildIntersectMatrix(IntersectMetric metric) {
        buildAllIntersectionsLazy()

        def matrix = new double[samples.length][samples.length]

        for (int i = 0; i < samples.length; i++) {
            matrix[i][i] = -1 // mask

            for (int j = i + 1; j < samples.length; j++) {
                matrix[i][j] = metric.value(pairedIntersectionsMat[i][j])
                matrix[j][i] = matrix[i][j]
            }
        }

        matrix
    }

    void buildAllIntersectionsLazy() {
        if (!allIntersectionsBuilt) {
            for (int i = 0; i < samples.length - 1; i++) {
                pairedIntersectionsMat[i][i + 1] = pairedIntersections[i]

                for (int j = i + 2; j < samples.length; j++) {
                    pairedIntersectionsMat[i][j] = intersectionUtil.generatePairedIntersection(samples[i], samples[j])
                }
            }
            allIntersectionsBuilt = true
        }
    }

    TimeCourse asTimeCourse() {
        def clonotypeMap = new HashMap<String, Clonotype[]>()

        pairedIntersections.eachWithIndex { PairedIntersection intersection, int sampleId ->
            intersection.clonotypes12.eachWithIndex { Clonotype clonotype12, int i ->
                def clonotype21 = intersection.clonotypes21[i]

                def key = intersectionUtil.generateKey(clonotype21)
                def entry = clonotypeMap[key]
                if (!entry)
                    clonotypeMap.put(key, entry = new Clonotype[samples.size()])
                entry[sampleId] = clonotype12
                entry[sampleId + 1] = clonotype21
            }
        }

        new TimeCourse(samples, clonotypeMap.values().collect { new DynamicClonotype(it) })
    }

    static final String HEADER = PairedIntersection.HEADER

    /*
    private double[] cumulative(String prop) {
        def d = 0
        def res = [0, pairedIntersections.collect {
            def dd = 1 - it."$prop"
            d += dd
        }].flatten()

        res.collect { it / d } as double[]
    }

    double[] getCumulativeF() {
        cumulative("f")
    }

    double[] getCumulativeD() {
        cumulative("d")
    }

    double[] getCumulativeR() {
        cumulative("r")
    } */


    @Override
    String toString() {
        pairedIntersections.join("\n")
    }
}
