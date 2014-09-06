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

class Spectratype {
    private final boolean aminoAcid, unweighted
    private final int min, max, len

    //private final static int SPECTRA_MIN = 8, SPECTRA_MAX = 32, SPECTRA_LEN = SPECTRA_MAX - SPECTRA_MIN + 1

    //private static int spectraBin(Clonotype clonotype) {
    //    Math.min(SPECTRA_MAX, Math.max(SPECTRA_MIN, clonotype.cdr3aa.length())) - SPECTRA_MIN
    //}
    private final double[] spectratype
    final String HEADER
    final int[] lengths

    Spectratype(boolean aminoAcid, boolean unweighted) {
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

    int bin(Clonotype clonotype) {
        Math.min(max, Math.max(min, (aminoAcid ? getCdr3AaLen(clonotype) : clonotype.cdr3nt.length()))) - min
    }

    private static int getCdr3AaLen(Clonotype clonotype) {
        (int) (clonotype.inFrame ? clonotype.cdr3aa.length() : (clonotype.cdr3nt.length() / 3))
    }

    public void addAll(Iterable<Clonotype> sample) {
        sample.each { Clonotype clonotype ->
            if (unweighted)
                this.spectratype[bin(clonotype)]++
            else
                this.spectratype[bin(clonotype)] += clonotype.freq
        }
    }

    public void clear() {
        for (int i = 0; i < len; i++)
            this.spectratype[i] = 0
    }

    double[] getHistogram() {
        double s = this.spectratype.collect().sum()

        def spectratype = new double[len]

        for (int i = 0; i < len; i++)
            spectratype[i] = this.spectratype[i] / s

        spectratype
    }

    @Override
    String toString() {
        histogram.collect().join("\t")
    }
}
