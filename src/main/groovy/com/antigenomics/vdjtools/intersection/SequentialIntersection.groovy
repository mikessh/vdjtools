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
    private final PairedIntersection[] intersectionSequence
    private final IntersectionUtil intersectionUtil

    /**
     * Builds a sequential intersection, i.e. a set of paired intersection of i-th and i+1-th samples
     * @param samples samples to intersect
     * @param intersectionUtil object used to build intersections
     */
    SequentialIntersection(Sample[] samples, IntersectionUtil intersectionUtil) {
        if (samples.size() < 3)
            throw new IllegalArgumentException("More than 2 samples should be provided")

        this.samples = samples
        this.intersectionUtil = intersectionUtil

        this.intersectionSequence = (0..(samples.length - 2)).collect { int i ->
            intersectionUtil.generatePairedIntersection(samples[i], samples[i + 1])
        } as PairedIntersection[]
    }

    /**
     * Creates sequential intersection from pairwise intersection matrix by cropping the +1 diagonal
     * @param pairedIntersectionMatrix base paired intersection matrix
     */
    SequentialIntersection(PairedIntersectionMatrix pairedIntersectionMatrix) {
        int n = pairedIntersectionMatrix.size() - 1
        this.samples = new Sample[n + 1]
        this.intersectionSequence = new PairedIntersection[n]
        this.intersectionUtil = pairedIntersectionMatrix.intersectionUtil

        for (int i = 0; i < n; i++) {
            intersectionSequence[i] = pairedIntersectionMatrix[i, i + 1]
            samples[i] = intersectionSequence[i].sample1
        }

        // last one
        samples[n] = intersectionSequence[n - 1].sample2
    }

    /**
     * Generates a time course for a given sequential intersection.
     * Only clonotypes met in at least two sequential samples are retained
     * @return clonotype abundance time course
     */
    TimeCourse asTimeCourse() {
        def clonotypeMap = new HashMap<String, Clonotype[]>()

        intersectionSequence.eachWithIndex { PairedIntersection intersection, int sampleIndex ->
            intersection.clonotypes12.eachWithIndex { Clonotype clonotype12, int i ->
                def clonotype21 = intersection.clonotypes21[i]

                def key = intersectionUtil.generateKey(clonotype21)
                def entry = clonotypeMap[key]
                if (!entry)
                    clonotypeMap.put(key, entry = new Clonotype[samples.size()])

                entry[sampleIndex] = clonotype12
                entry[sampleIndex + 1] = clonotype21
            }
        }

        new TimeCourse(samples, clonotypeMap.values().collect { new DynamicClonotype(it) })
    }

    @Override
    String toString() {
        "Samples=" + samples.join(",")
    }

    void print(PrintWriter pw, boolean includeHeader) {
        if (includeHeader)
            pw.println("#" + PairedIntersection.HEADER)

        pw.println(intersectionSequence.join("\n"))
    }
}
