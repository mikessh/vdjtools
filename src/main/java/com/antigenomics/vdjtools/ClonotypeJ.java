/**
 * Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.vdjtools;

import java.util.HashSet;
import java.util.Set;

public class ClonotypeJ {
    private final int count;
    private final double freq;

    private final int[] segmPoints;
    private final String v, d, j;
    private final String cdr1nt, cdr2nt, cdr3nt, cdr1aa, cdr2aa, cdr3aa;

    private final boolean inFrame, isComplete, noStop;

    private final Set<Mutation> mutations;

    public ClonotypeJ(int count, double freq,
                      int[] segmPoints, String v, String d, String j,
                      String cdr1nt, String cdr2nt, String cdr3nt,
                      String cdr1aa, String cdr2aa, String cdr3aa,
                      boolean inFrame, boolean isComplete, boolean noStop,
                      Set<Mutation> mutations) {
        this.count = count;
        this.freq = freq;
        this.segmPoints = segmPoints;
        this.v = v;
        this.d = d;
        this.j = j;
        this.cdr1nt = cdr1nt;
        this.cdr2nt = cdr2nt;
        this.cdr3nt = cdr3nt;
        this.cdr1aa = cdr1aa;
        this.cdr2aa = cdr2aa;
        this.cdr3aa = cdr3aa;
        this.inFrame = inFrame;
        this.isComplete = isComplete;
        this.noStop = noStop;
        this.mutations = mutations;
    }

    public int getCount() {
        return count;
    }

    public double getFreq() {
        return freq;
    }

    public String getV() {
        return v;
    }

    public String getD() {
        return d;
    }

    public String getJ() {
        return j;
    }

    public String getCdr1nt() {
        return cdr1nt;
    }

    public String getCdr3nt() {
        return cdr3nt;
    }

    public String getCdr2nt() {
        return cdr2nt;
    }

    public String getCdr1aa() {
        return cdr1aa;
    }

    public String getCdr2aa() {
        return cdr2aa;
    }

    public String getCdr3aa() {
        return cdr3aa;
    }

    public boolean inFrame() {
        return inFrame;
    }

    public boolean isComplete() {
        return isComplete;
    }

    public boolean noStop() {
        return noStop;
    }

    public int getVEnd() {
        return segmPoints[0];
    }

    public int getDStart() {
        return segmPoints[1];
    }

    public int getDEnd() {
        return segmPoints[2];
    }

    public int getJStart() {
        return segmPoints[3];
    }

    public int getVDIns() {
        return segmPoints[1] >= 0 ? segmPoints[1] - segmPoints[0] + 1 : -1;
    }

    public int getDJIns() {
        return segmPoints[2] >= 0 ? segmPoints[3] - segmPoints[2] + 1 : -1;
    }

    public int getVJIns() {
        return segmPoints[1] >= 0 && segmPoints[2] >= 0 ?
                getVDIns() + getDJIns() :
                segmPoints[3] - segmPoints[1] + 1;
    }

    public String getKey() {
        StringBuilder key = new StringBuilder(v + "_" + cdr3nt + "_" + j + "_");
        for (Mutation mutation : mutations) {
            key.append("_").append(mutation)
        }
    }
}
