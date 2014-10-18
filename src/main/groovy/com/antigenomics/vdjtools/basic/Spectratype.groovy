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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.intersection.IntersectionType

class Spectratype {
    private final boolean aminoAcid, unweighted
    private final int min, max, len

    private final double[] spectratype
    final String HEADER
    final int[] lengths

    private int count = 0
    private double freq = 0

    public Spectratype(IntersectionType intersectionType, boolean unweighted) {
        this(intersectionType.aminoAcid, unweighted)
    }

    public Spectratype(boolean aminoAcid, boolean unweighted) {
        this.aminoAcid = aminoAcid
        this.unweighted = unweighted

        if (aminoAcid) {
            this.min = 7
            this.max = 23
        } else {
            this.min = 21
            this.max = 69
        }

        this.len = max - min + 1
        this.spectratype = new double[len]
        this.lengths = (min..max) as double[]
        this.HEADER = lengths.collect().join("\t")
    }

    private int bin(Clonotype clonotype) {
        Math.min(max, Math.max(min, (aminoAcid ? getCdr3AaLen(clonotype) : clonotype.cdr3nt.length()))) - min
    }

    private static int getCdr3AaLen(Clonotype clonotype) {
        (int) (clonotype.inFrame ? clonotype.cdr3aa.length() : (clonotype.cdr3nt.length() / 3))
    }

    public void add(Clonotype clonotype) {
        if (unweighted) {
            this.spectratype[bin(clonotype)]++
            freq++
        } else {
            this.spectratype[bin(clonotype)] += clonotype.freq
            freq += clonotype.freq
        }
        count++
    }

    public void addAll(Iterable<Clonotype> sample) {
        sample.each { Clonotype clonotype ->
            add(clonotype)
        }
    }

    public void addAll(Spectratype other) {
        if (other.aminoAcid != this.aminoAcid || other.unweighted != this.unweighted)
            throw new Exception("Can't add Spectratype of different type")

        this.count += other.count
        this.freq += other.freq
        other.histogram.eachWithIndex { it, ind ->
            histogram[ind] += it
        }
    }

    public List<Clonotype> addAllFancy(Iterable<Clonotype> sample, int top) {
        def topClonotypes = new LinkedList<Clonotype>();
        int counter = 0
        sample.each { Clonotype clonotype ->
            if (counter++ < top)
                topClonotypes.add(clonotype)
            else
                add(clonotype)
        }
        topClonotypes.sort { -it.freq }
    }

    public void clear() {
        for (int i = 0; i < len; i++)
            this.spectratype[i] = 0
        count = 0
        freq = 0
    }

    public double[] getHistogram() {
        getHistogram(true)
    }

    public double[] getHistogram(boolean normalized) {
        def spectratype = new double[len]
        def _freq = normalized ? freq : 1.0

        for (int i = 0; i < len; i++)
            spectratype[i] = this.spectratype[i] / _freq

        spectratype
    }

    public int getLen() {
        return len
    }

    public double getFreq() {
        return freq
    }

    public int getCount() {
        return count
    }

    @Override
    String toString() {
        histogram.collect().join("\t")
    }
}
