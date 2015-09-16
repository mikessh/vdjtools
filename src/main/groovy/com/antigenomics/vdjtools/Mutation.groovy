/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.sample.Clonotype

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

    //static Mutation cdr3Mutation(int ntPos, char fromNt, char toNt,
    //                             int aaPos, char fromAa, char toAa,
    //                             boolean directed, Clonotype parent) {
    //    new Mutation("CDR3",
    //            ntPos, aaPos,
    //            fromAa, toAa, fromNt, toNt,
    //            directed, parent)
    //}

    String getV() {
        parent.v
    }

    double getFreq() {
        directed ? parent.freq : altFreq
    }

    //String getMotif(int lSize, int rSize) {
    //    int from = ntPos - lSize, to = ntPos + rSize + 1
    //    region == "CDR3" ? parent.getSubSequence(from, to) :
    //            parent.getVSegmentData().getSubSequence(from, to)
    //}

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
