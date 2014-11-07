/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 *
 * Last modified on 30.10.2014 by mikesh
 */

package com.antigenomics.vdjtools;

import com.antigenomics.vdjtools.sample.Sample;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class Clonotype implements Comparable<Clonotype>, Countable {
    private final Sample parent;
    private final int count;
    private final double freq;
    private final String key;

    private final int[] segmPoints;
    private final String v, d, j;
    private final String cdr1nt, cdr2nt, cdr3nt,
            cdr1aa, cdr2aa, cdr3aa;

    private final boolean inFrame, isComplete, noStop;

    private final Set<Mutation> mutations;

    public Clonotype(Sample parent, int count, double freq,
                     int[] segmPoints, String v, String d, String j,
                     String cdr1nt, String cdr2nt, String cdr3nt,
                     String cdr1aa, String cdr2aa, String cdr3aa,
                     boolean inFrame, boolean isComplete, boolean noStop,
                     Set<Mutation> mutations) {
        this.parent = parent;
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

        StringBuilder key = new StringBuilder(v).append(KEY_SEP).
                append(cdr3nt).append(KEY_SEP).
                append(j).append(KEY_SEP);

        if (mutations != null)
            for (Mutation mutation : mutations) {
                key.append(mutation.getNtString()).append(MUT_SEP);
            }

        key.setLength(key.length() - 1);

        this.key = key.toString();
    }

    public Clonotype(Clonotype toClone) {
        this(toClone, toClone.parent, toClone.count);
    }

    public Clonotype(Clonotype toClone, Sample newParent) {
        this(toClone, newParent, toClone.count);
    }

    public Clonotype(Clonotype toClone, Sample newParent, int newCount) {
        this(newParent, newCount, toClone.freq,
                toClone.segmPoints, toClone.v, toClone.d, toClone.j,
                toClone.cdr1nt, toClone.cdr2nt, toClone.cdr3nt,
                toClone.cdr1aa, toClone.cdr2aa, toClone.cdr3aa,
                toClone.inFrame, toClone.isComplete, toClone.noStop,
                toClone.mutations != null ? new HashSet<Mutation>() : null);

        if (toClone.mutations != null)
            for (Mutation mutation : toClone.mutations)
                mutations.add(mutation.reassignParent(this));
    }

    public int getCount() {
        return count;
    }

    public double getFreqAsInInput() {
        return freq;
    }

    public double getFreq() {
        return count / (double) parent.getCount();
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

    // todo: extract
    public String getCdr1nt() {
        return cdr1nt;
    }

    public String getCdr2nt() {
        return cdr2nt;
    }

    public String getCdr3nt() {
        return cdr3nt;
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

    public boolean isInFrame() {
        return inFrame;
    }

    public boolean isComplete() {
        return isComplete;
    }

    public boolean isNoStop() {
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

    public int getInsertSize() {
        return segmPoints[1] >= 0 && segmPoints[2] >= 0 ? getVDIns() + getDJIns() : -1;
    }

    public int getNDNSize() {
        return segmPoints[3] - segmPoints[1] + 1;
    }

    public String getBlank() {
        return ".";
    }

    public Sample getParent() {
        return parent;
    }

    public Set<Mutation> getMutations() {
        return mutations != null ? Collections.unmodifiableSet(mutations) : null;
    }

    private static final String KEY_SEP = "_", MUT_SEP = "|";

    public String getKey() {
        return key;
    }

    @Override
    public int compareTo(Clonotype o) {
        return -Integer.compare(this.count, o.count);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Clonotype that = (Clonotype) o;

        return key.equals(that.key) && parent.equals(that.parent);
    }

    @Override
    public int hashCode() {
        return 31 * parent.hashCode() + key.hashCode();
    }
}
