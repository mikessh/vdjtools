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
import com.antigenomics.vdjtools.sample.SampleUtil
import com.antigenomics.vdjtools.timecourse.DynamicClonotype
import com.antigenomics.vdjtools.timecourse.TimeCourse

class PairedIntersection {
    final Sample sample1, sample2

    final int uniq12, count12, count21
    final double freq12, freq21

    private final List<Clonotype> clonotypes12, clonotypes21
    private Double r = null

    PairedIntersection(Sample sample1, Sample sample2,
                             int uniq12,
                             int count12, int count21,
                             double freq12, double freq21,
                             List<Clonotype> clonotypes12, List<Clonotype> clonotypes21) {
        this.sample1 = sample1
        this.sample2 = sample2
        this.uniq12 = uniq12
        this.count12 = count12
        this.count21 = count21
        this.freq12 = freq12
        this.freq21 = freq21
        this.clonotypes12 = clonotypes12
        this.clonotypes21 = clonotypes21
    }

    /**
     * Gets the correlation between log-frequencies of clonotypes, detected in both samples
     * @return
     */
    double getCorrelation() {
        // lazy compute
        r ?: (r = SampleUtil.correlation(clonotypes12, clonotypes21))
    }

    TimeCourse asTimeCourse() {
        def dynamicClontypes = new ArrayList<DynamicClonotype>()
        for (int i = 0; i < uniq12; i++)
            dynamicClontypes.add(new DynamicClonotype([clonotypes12[i], clonotypes21[i]] as Clonotype[]))
        def samples = [sample1, sample2] as Sample[]
        new TimeCourse(samples, dynamicClontypes)
    }

    List<Clonotype> getClonotypes12() {
        return clonotypes12
    }

    List<Clonotype> getClonotypes21() {
        return clonotypes21
    }

    //
    // Print
    //

    final static String HEADER = "clones1\tclones2\tclones12\t" +
            "count1\tcount2\tcount12\tcount21\t" +
            "freq1\tfreq2\tfreq12\tfreq21\t" +
            "R"

    @Override
    String toString() {
        [sample1.diversity, sample2.diversity, uniq12,
         sample1.count, sample2.count, count12, count21,
         sample1.freq, sample2.freq, freq12, freq21,
         getCorrelation()].join("\t")
    }
}
