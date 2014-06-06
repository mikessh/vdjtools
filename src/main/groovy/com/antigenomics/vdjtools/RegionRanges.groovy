package com.antigenomics.vdjtools

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
class RegionRanges {
    final Range fr1, cdr1, fr2, cdr2, fr3

    RegionRanges(String fr1from, String fr1to,
                 String cdr1from, String cdr1to,
                 String fr2from, String fr2to,
                 String cdr2from, String cdr2to,
                 String fr3from, String fr3to) {
        fr1 = new Range(fr1from, fr1to)
        cdr1 = new Range(cdr1from, cdr1to)
        fr2 = new Range(fr2from, fr2to)
        cdr2 = new Range(cdr2from, cdr2to)
        fr3 = new Range(fr3from, fr3to)
    }

    int posToRegion(int pos) {
        if (fr1.inRange(pos))
            return 0
        if (cdr1.inRange(pos))
            return 1
        if (fr2.inRange(pos))
            return 2
        if (cdr2.inRange(pos))
            return 3
        if (fr3.inRange(pos))
            return 4
        -1
    }

    static String regionId2Name(int regionId) {
        switch (regionId) {
            case 0:
                return "FR1"
            case 1:
                return "CDR1"
            case 2:
                return "FR2"
            case 3:
                return "CDR2"
            case 4:
                return "FR3"
            default:
                return "N/A"
        }
    }

    static final String HEADER = (0..<N_REGIONS).collect { regionId2Name(it) }.join("\t")
    static final int N_REGIONS = 6
}
