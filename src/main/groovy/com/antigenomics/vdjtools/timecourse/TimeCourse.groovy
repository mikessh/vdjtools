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

package com.antigenomics.vdjtools.timecourse

import com.antigenomics.vdjtools.sample.Sample

class TimeCourse implements Iterable<DynamicClonotype> {
    final Sample[] samples
    private final List<DynamicClonotype> clonotypes

    private boolean sorted = false

    TimeCourse(Sample[] samples, Collection<DynamicClonotype> clonotypes) {
        this.samples = samples
        this.clonotypes = clonotypes
    }

    CollapsedTimeCourse collapseBelow(int top) {
        lazySort()

        final double[] upperFrequency = new double[samples.length],
                       lowerFrequency = new double[samples.length],
                       remainingFrequency = new double[samples.length]

        (0..<top).each { int index ->
            clonotypes[index].frequencies.eachWithIndex { double f, int i ->
                upperFrequency[i] += f
            }
        }

        (top..<clonotypes.size()).each { int index ->
            clonotypes[index].frequencies.eachWithIndex { double f, int i ->
                lowerFrequency[i] += f
            }
        }

        for (int i = 0; i < samples.size(); i++)
            remainingFrequency[i] = samples[i].freq - upperFrequency[i] - lowerFrequency[i]

        new CollapsedTimeCourse(samples, clonotypes[0..<top], upperFrequency, lowerFrequency, remainingFrequency)
    }

    private void lazySort() {
        if (!sorted) {
            clonotypes.sort { -it.meanFrequency }
            sorted = true
        }
    }

    String getHeader() {
        [DynamicClonotype.PRINT_FIELDS, samples.collect { it.metadata.sampleId }].flatten().join("\t")
    }

    void print(PrintWriter pw, boolean addHeader) {
        lazySort()

        if (addHeader)
            pw.println(header)

        clonotypes.each {
            pw.println(it.toString() + "\t" + it.frequencies.collect().join("\t"))
        }
    }

    @Override
    Iterator<DynamicClonotype> iterator() {
        clonotypes.iterator()
    }
}
