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

class SpectratypeV {
    private final Map<String, Spectratype> spectratypes = new HashMap<>()
    private final boolean aminoAcid, unweighted

    SpectratypeV(IntersectionType intersectionType, boolean unweighted) {
        this(intersectionType.aminoAcid, unweighted)
    }

    public SpectratypeV(boolean aminoAcid, boolean unweighted) {
        this.aminoAcid = aminoAcid
        this.unweighted = unweighted
    }

    public void clear() {
        spectratypes.values().each { it.clear() }
    }

    public void addAll(Iterable<Clonotype> sample) {
        sample.each { clonotype ->
            def vSpectra = spectratypes[clonotype.v]
            if (!vSpectra)
                spectratypes.put(clonotype.v, vSpectra = new Spectratype(aminoAcid, unweighted))
            vSpectra.add(clonotype)
        }
    }

    public List<String> vSegmentList() {
        spectratypes.sort { -it.value.freq }.collect { it.key }
    }

    public Spectratype getAt(String vSegmentName) {
        spectratypes[vSegmentName]
    }

    public Map<String, Spectratype> collapse(int top) {
        def collapsedSpectratypes = new HashMap(), otherSpectratype
        collapsedSpectratypes.put("other", otherSpectratype = new Spectratype(aminoAcid, unweighted))

        spectratypes.sort { -it.value.freq }.eachWithIndex { it, ind ->
            if (ind < top)
                collapsedSpectratypes.put(it.key, it.value)
            else
                otherSpectratype.addAll(it.value)
        }

        collapsedSpectratypes
    }
}
