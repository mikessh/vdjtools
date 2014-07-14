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

package com.antigenomics.vdjtools.spectratype

import com.antigenomics.vdjtools.Clonotype

class SpectratypePeak {
    final String v, signature
    final int cdr3length
    final List<Clonotype> clonotypes = new ArrayList<>()
    int clones = 0
    double freq = 0.0

    SpectratypePeak(Clonotype clonotype) {
        this.v = clonotype.v
        this.cdr3length = clonotype.cdr3aa.length()
        this.signature = v + "\t" + cdr3length
    }

    void append(Clonotype clonotype) {
        clones++
        freq += clonotype.freq
        clonotypes.add(clonotype)
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        SpectratypePeak that = (SpectratypePeak) o

        if (signature != that.signature) return false

        true
    }

    @Override
    int hashCode() {
        signature.hashCode()
    }

    final static String HEADER = "v\tlen\tuniq\tfreq"

    @Override
    String toString() {
        [signature, clones, freq].join("\t")
    }
}
