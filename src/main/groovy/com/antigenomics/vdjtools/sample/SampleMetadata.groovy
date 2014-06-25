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

class SampleMetadata {
    final String fileName, sampleId
    final List<String> entries

    SampleMetadata(String fileName, String sampleId, List<String> entries) {
        this.fileName = fileName
        this.sampleId = sampleId
        this.entries = entries
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        SampleMetadata that = (SampleMetadata) o

        if (sampleId != that.sampleId) return false

        return true
    }

    int hashCode() {
        return sampleId.hashCode()
    }
}
