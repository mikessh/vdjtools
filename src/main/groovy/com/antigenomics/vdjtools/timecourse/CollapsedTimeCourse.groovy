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

class CollapsedTimeCourse extends TimeCourse {
    //final Sample[] samples
    //final List<DynamicClonotype> clonotypes
    final double[] upperFrequency, lowerFrequency, remainingFrequency

    CollapsedTimeCourse(Sample[] samples, List<DynamicClonotype> clonotypes,
                        double[] upperFrequency, double[] lowerFrequency, double[] remainingFrequency) {
        super(samples, clonotypes)
        this.upperFrequency = upperFrequency
        this.lowerFrequency = lowerFrequency
        this.remainingFrequency = remainingFrequency
    }


    @Override
    void print(PrintWriter pw, boolean addHeader) {
        super.print(pw, addHeader)

        pw.println(DynamicClonotype.PRINT_FIELDS.collect {
            "Not-shown"
        }.join("\t") + "\t-1\t" + lowerFrequency.collect().join("\t"))


        pw.println(DynamicClonotype.PRINT_FIELDS.collect {
            "Non-overlapping"
        }.join("\t") + "\t-2\t" + remainingFrequency.collect().join("\t"))
    }
}
