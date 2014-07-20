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

package com.antigenomics.vdjtools.sample

class SamplePair {
    final Sample sample1, sample2

    SamplePair(Sample sample1, Sample sample2) {
        this.sample1 = sample1
        this.sample2 = sample2
    }

    Sample getAt(int index) {
        switch (index) {
            case 0:
                return sample1
            case 1:
                return sample2
        }
        throw new IllegalArgumentException("Index must be either 0 or 1")
    }
}
