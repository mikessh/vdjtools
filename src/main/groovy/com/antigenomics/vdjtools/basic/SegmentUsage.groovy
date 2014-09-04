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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.Clonotype

class SegmentUsage {
    private final Map<String, double[]> vSegmentUsage, jSegmentUsage
    private final int n
    private final boolean unweighted

    public SegmentUsage(int nSamples, boolean unweighted) {
        this.n = nSamples
        this.unweighted = unweighted
    }

    public void addAll(Iterable<Clonotype> sample, int sampleIndex) {
        sample.each { Clonotype clonotype ->
            def vArray = vSegmentUsage[clonotype.v],
                jArray = jSegmentUsage[clonotype.j]

            if (!vArray)
                vSegmentUsage.put(clonotype.v, vArray = new double[n])

            if (!jArray)
                jSegmentUsage.put(clonotype.j, jArray = new double[n])

            def increment = unweighted ? 1 : clonotype.freq

            vArray[sampleIndex] += increment
            jArray[sampleIndex] += increment
        }
    }
}
