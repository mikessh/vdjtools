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

class IntersectionResult {
    final int clones1, clones2, intersectionClones
    final double freq1, freq2, intersectionFreq

    IntersectionResult(int clones1, int clones2, int intersectionClones,
                     double freq1, double freq2, double intersectionFreq) {
        this.clones1 = clones1
        this.clones2 = clones2
        this.intersectionClones = intersectionClones
        this.freq1 = freq1
        this.freq2 = freq2
        this.intersectionFreq = intersectionFreq
    }
}
