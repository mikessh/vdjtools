package com.antigenomics.vdjtools.segment

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

class SegmentData {
    final String segmentName
    final String segmentSequence
    final List<Range> regionMarkup

    SegmentData(String segmentName, String segmentSequence, List<Range> regionMarkup) {
        this.segmentName = segmentName
        this.segmentSequence = segmentSequence
        this.regionMarkup = regionMarkup
    }

    int deduceRegion(int pos) {
        regionMarkup.findIndexOf { it.contains(pos)}
    }
}
