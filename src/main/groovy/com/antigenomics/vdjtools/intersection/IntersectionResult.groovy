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

class IntersectionResult {
    final int uniq1, uniq2, uniq12
    final Collection<Clonotype> clonotypes1, clonotypes2
    final double freq1, freq2, freq12

    IntersectionResult(int uniq1, int uniq2, int uniq12,
                       double freq1, double freq2, double freq12,
                       Collection<Clonotype> clonotypes1, Collection<Clonotype> clonotypes2) {
        this.uniq1 = uniq1
        this.uniq2 = uniq2
        this.uniq12 = uniq12
        this.freq1 = freq1
        this.freq2 = freq2
        this.freq12 = freq12
        this.clonotypes1 = clonotypes1
        this.clonotypes2 = clonotypes2
    }

    final static String HEADER = "uniq1\tuniq2\tuniq12\tfreq1\tfreq2\tfreq12"

    @Override
    String toString() {
        [uniq1, uniq2, uniq12, freq1, freq2, freq12].join("\t")
    }
}
