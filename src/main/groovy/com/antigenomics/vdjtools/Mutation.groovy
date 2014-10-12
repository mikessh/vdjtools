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

package com.antigenomics.vdjtools

class Mutation {
    final String key, ntString, aaString, region
    final int ntPos, aaPos
    final char fromAa, toAa, fromNt, toNt
    final boolean isSilent, isStop, directed
    final Clonotype parent
    double altFreq

    public Mutation(String region,
                     int ntPos, int aaPos,
                     char fromAa, char toAa, char fromNt, char toNt,
                     boolean directed,
                     Clonotype parent) {
        this.parent = parent
        this.directed = directed
        def delim = directed ? ">" : "<>"
        this.isStop = toAa == "*" || (!directed && fromAa == "*")
        this.isSilent = fromAa == toAa
        this.toNt = toNt
        this.fromNt = fromNt
        this.toAa = toAa
        this.fromAa = fromAa
        this.aaPos = aaPos
        this.ntPos = ntPos
        this.region = region
        this.aaString = aaPos + delim + fromAa + ":" + toAa
        this.ntString = ntPos + delim + fromNt + ":" + toNt
        this.key = region + ":" + ntString
    }

    Mutation reassignParent(Clonotype parent) {
        new Mutation(region,
                ntPos, aaPos,
                fromAa, toAa, fromNt, toNt,
                directed, parent)
    }

    static Mutation cdr3Mutation(int ntPos, char fromNt, char toNt,
                                 int aaPos, char fromAa, char toAa,
                                 boolean directed, Clonotype parent) {
        new Mutation("CDR3",
                ntPos, aaPos,
                fromAa, toAa, fromNt, toNt,
                directed, parent)
    }

    String getV() {
        parent.v
    }

    double getFreq() {
        directed ? parent.freq : altFreq
    }

    String getMotif(int lSize, int rSize) {
        int from = ntPos - lSize, to = ntPos + rSize + 1
        region == "CDR3" ? parent.getSubSequence(from, to) :
                parent.getVSegmentData().getSubSequence(from, to)
    }

    //
    // Intersection
    //

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Mutation that = (Mutation) o

        if (key != that.key) return false

        return true
    }

    @Override
    int hashCode() {
        return key.hashCode()
    }

    //
    // Display
    //

    final static String HEADER = "display_name\t" +
            "nt_shm\taa_shm\tregion\tsilent\tdirected"

    String getDisplayName() {
        isSilent ? "S" : (region + ":" + fromAa + ">" + toAa)
    }

    @Override
    String toString() {
        [displayName, ntString, aaString, region, isSilent, directed].join("\t")
    }
}
