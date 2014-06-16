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
class SegmentUtil {
    static String regionId2Name(int regionId) {
        switch (regionId) {
            case 0:
                return "FW1"
            case 1:
                return "CDR1"
            case 2:
                return "FW2"
            case 3:
                return "CDR2"
            case 4:
                return "FW3"
            case 5:
                return "CDR3"
            default:
                return "N/A"
        }
    }

    static int regionName2Id(String regionName) {
        switch (regionName) {
            case "FW1":
                return 0
            case "CDR1":
                return 1
            case "FW2":
                return 2
            case "CDR2":
                return 3
            case "FW3":
                return 4
            case "CDR3":
                return 5
            default:
                return -1
        }
    }

    static final int N_REGIONS = 6
}
