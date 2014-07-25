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
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.timecourse.DynamicClonotype
import com.antigenomics.vdjtools.timecourse.TimeCourse

class PairedIntersection {
    final SamplePair parent

    final int div12, count12, count21
    final double freq12, freq21, freq12e, freq21e, freq12p, freq21p

    private final List<Clonotype> clonotypes12, clonotypes21
    private final double r

    PairedIntersection(SamplePair parent,
                       int div12, int count12, int count21,
                       double freq12, double freq21,
                       double freq12e, double freq21e,
                       double freq12p, double freq21p,
                       double r,
                       List<Clonotype> clonotypes12, List<Clonotype> clonotypes21) {
        this.parent = parent
        this.div12 = div12
        this.count12 = count12
        this.count21 = count21
        this.freq12 = freq12
        this.freq21 = freq21
        this.freq12e = freq12e
        this.freq21e = freq21e
        this.freq12p = freq12p
        this.freq21p = freq21p
        this.r = r
        this.clonotypes12 = clonotypes12
        this.clonotypes21 = clonotypes21
    }

    /**
     * Gets the correlation between log-frequencies of clonotypes that were detected in both samples
     * @return correlation
     */
    double getR() {
        r
        // lazy compute
        //r ?: (r = SampleUtil.correlation(clonotypes12, clonotypes21))
    }

    double getNormR() {
        (1 + r) / 2
    }

    /**
     * Gets the geometric mean of intersection frequencies
     * @return normalized frequency
     */
    double getF() {
        Math.sqrt(freq12 * freq21)
    }

    double getNormF() {
        getF()
    }

    /**
     * Gets the diversity of intersection, normalized to its expected value
     * @return normalized diversity
     */
    double getD() {
        div12 / Math.sqrt((double) div1 * (double) div2)
    }

    double getNormD() {
        def maxValue = div1 < div2 ? Math.sqrt(div1 / div2) : Math.sqrt(div2 / div1)
        getD() / maxValue
    }

    TimeCourse asTimeCourse() {
        def dynamicClontypes = new ArrayList<DynamicClonotype>()
        for (int i = 0; i < div12; i++)
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

    Sample getSample1() {
        parent.sample1
    }

    Sample getSample2() {
        parent.sample2
    }

    int getDiv1() {
        sample1.div
    }

    int getDiv2() {
        sample2.div
    }

    int getCount1() {
        sample1.count
    }

    int getCount2() {
        sample2.count
    }

    double getFreq1() {
        sample1.freq
    }

    double getFreq2() {
        sample2.freq
    }

    //
    // Print
    //

    final static String HEADER = "1_sample_id\t2_sample_id\t" +
            "div1\tdiv2\tdiv12\t" +
            "count1\tcount2\tcount12\tcount21\t" +
            "freq1\tfreq2\tfreq12\tfreq21\t" +
            "freq12e\tfreq21e\tfreq12p\tfreq21p\t" +
            "F\tD\tR"

    @Override
    String toString() {
        [sample1.metadata.sampleId, sample2.metadata.sampleId,
         div1, div2, div12,
         count1, count2, count12, count21,
         freq1, freq2, freq12, freq21,
         freq12e, freq21e, freq12p, freq21p,
         getNormF(), getNormD(), getNormR()].join("\t")
    }
}
