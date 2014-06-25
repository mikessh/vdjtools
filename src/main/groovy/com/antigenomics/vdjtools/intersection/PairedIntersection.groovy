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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.Clonotype

class PairedIntersection {
    final IntersectionType intersectionType
    private final Collection<ClonotypeWrapper> sample1, sample2

    PairedIntersection(Collection<Clonotype> sample1, Collection<Clonotype> sample2,
                       IntersectionType intersectionType) {
        this.intersectionType = intersectionType
        this.sample1 = sample1.collect { new ClonotypeWrapper(it, intersectionType) }
        this.sample2 = sample2.collect { new ClonotypeWrapper(it, intersectionType) }
    }

    IntersectionResult intersect() {
        def intersection = new HashMap<ClonotypeWrapper, ClonotypeWrapper>()

        boolean flip = sample1.size() < sample2.size()

        Collection<ClonotypeWrapper> _sample1, _sample2

        (_sample1, _sample2) = flip ? [sample2, sample1] : [sample1, sample2]

        double freq1 = 0, freq2 = 0, intersectedFreq = 0
        int clones1 = _sample1.size(), clones2 = _sample2.size(), intersectedClones = 0

        _sample1.each {
            intersection.put(it, it)
            freq1 += it.clonotype.freq
        }

        _sample2.each {
            def other = intersection[it]

            freq2 += it.clonotype.freq

            if (other != null) {
                intersectedClones++
                intersectedFreq += Math.sqrt(it.clonotype.freq * other.clonotype.freq)
            }
        }

        flip ? new IntersectionResult(clones2, clones1, intersectedClones, freq2, freq1, intersectedFreq) :
                new IntersectionResult(clones1, clones2, intersectedClones, freq1, freq2, intersectedFreq)
    }

    private class ClonotypeWrapper {
        final Clonotype clonotype
        final IntersectionType intersectionType

        ClonotypeWrapper(Clonotype clonotype,
                         IntersectionType intersectionType) {
            this.clonotype = clonotype
            this.intersectionType = intersectionType
        }

        boolean equals(o) {
            if (this.is(o)) return true

            ClonotypeWrapper clonotypeWrapper = (ClonotypeWrapper) o

            switch (intersectionType) {
                case IntersectionType.Nucleotide:
                    if (clonotype.cdr3nt != clonotypeWrapper.clonotype.cdr3nt)
                        return false
                    break
                case IntersectionType.AminoAcid:
                    if (clonotype.cdr3aa != clonotypeWrapper.clonotype.cdr3aa)
                        return false
                    break
                case IntersectionType.AminoAcidNonNucleotide:
                    if (clonotype.cdr3aa != clonotypeWrapper.clonotype.cdr3aa ||
                            clonotype.cdr3nt == clonotypeWrapper.clonotype.cdr3nt)
                        return false
                    break
            }

            return true
        }

        int hashCode() {
            intersectionType == IntersectionType.Nucleotide ?
                    clonotype.cdr3nt.hashCode() : clonotype.cdr3aa.hashCode()
        }
    }
}