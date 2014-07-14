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
import com.antigenomics.vdjtools.ClonotypeUtil

class IntersectionResult {
    final int uniq1, uniq2, uniq12
    final List<Clonotype> clonotypes12, clonotypes21
    final double freq1, freq2, freq12, freq21
    private Double r = null

    IntersectionResult(int uniq1, int uniq2, int uniq12,
                       double freq1, double freq2, double freq12, double freq21,
                       List<Clonotype> clonotypes12, List<Clonotype> clonotypes21) {
        this.uniq1 = uniq1
        this.uniq2 = uniq2
        this.uniq12 = uniq12
        this.freq1 = freq1
        this.freq2 = freq2
        this.freq12 = freq12
        this.freq21 = freq21
        this.clonotypes12 = clonotypes12
        this.clonotypes21 = clonotypes21
    }

    double getCorrelation() {
        r ?: (r = ClonotypeUtil.correlation(clonotypes12, clonotypes21))
    }

    double getMeanFrequency() {
        Math.sqrt(freq21 * freq12)
    }

    final static String HEADER = "uniq1\tuniq2\tuniq12\tfreq1\tfreq2\tfreq12\tfreq21\tR"

    @Override
    String toString() {
        [uniq1, uniq2, uniq12, freq1, freq2, freq12, freq21].join("\t")
    }
}
