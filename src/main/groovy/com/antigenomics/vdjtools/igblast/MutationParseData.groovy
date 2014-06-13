package com.antigenomics.vdjtools.igblast

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

class MutationParseData {
    final String key, aaString, region
    final ntPos, aaPos
    final char fromAA, toAA, fromNT, toNT
    final boolean isSilent, isStop
    final int nReads, nEvents
    final double fReads, fEvents

    MutationParseData(String mutationString) {
        def splitString = mutationString.split(",")

        def splitCountString = splitString[0].split(":")
        nReads = splitCountString[0].toInteger()
        fReads = splitCountString[1].toDouble()
        nEvents = splitCountString[2].toInteger()
        fEvents = splitCountString[3].toDouble()

        key = splitString[1]
        def splitNTString = key.split("[:>]")
        ntPos = splitNTString[0].toInteger()
        fromNT = splitNTString[1]
        toNT = splitNTString[2]

        aaString = splitString[2]
        def splitAAString = aaString.split("[:>]")
        aaPos = splitAAString[0].toInteger()
        fromAA = splitAAString[1]
        toAA = splitAAString[2]
        isSilent = fromAA == toAA
        isStop = toAA == "*"

        region = splitString[3]
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        MutationParseData that = (MutationParseData) o

        if (key != that.key) return false

        return true
    }

    int hashCode() {
        return key.hashCode()
    }

    //
    // Display
    //

    final static String EDGE_HEADER = "display_name\t" +
            "nt_shm\taa_shm\tregion\tsilent"

    String getDisplayName() {
        isSilent ? "S" : (region + ":" + fromAA + ">" + toAA)
    }

    @Override
    String toString() {
        [displayName, key, aaString, region, isSilent].join("\t")
    }
}
