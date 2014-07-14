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
    private final Collection<ClonotypeWrapper> sample1, sample2

    final IntersectionType intersectionType

    PairedIntersection(Collection<Clonotype> sample1, Collection<Clonotype> sample2,
                       IntersectionType intersectionType) {
        this.intersectionType = intersectionType
        this.sample1 = sample1.collect { new ClonotypeWrapper(it, intersectionType) }
        this.sample2 = sample2.collect { new ClonotypeWrapper(it, intersectionType) }
    }

    IntersectionResult intersect(boolean store) {
        def intersection = new HashMap<ClonotypeWrapper, ClonotypeWrapper>()

        boolean flip = sample1.size() < sample2.size()

        Collection<ClonotypeWrapper> _sample1, _sample2

        (_sample1, _sample2) = flip ? [sample2, sample1] : [sample1, sample2]

        double freq1 = 0, freq2 = 0, freq12 = 0, freq21 = 0
        int clones1 = _sample1.size(), clones2 = _sample2.size(), clones12 = 0

        _sample1.each {
            intersection.put(it, it)
            freq1 += it.clonotype.freq
        }

        final List<Clonotype> clonotypes1 = new ArrayList<>(),
                clonotypes2 = new ArrayList<>()

        _sample2.each {
            def other = intersection[it]

            freq2 += it.clonotype.freq

            if (other != null) {
                clones12++
                freq12 += it.clonotype.freq
                freq21 += other.clonotype.freq

                if (store) {
                    clonotypes1.add(other.clonotype)
                    clonotypes2.add(it.clonotype)
                }
            }
        }

        flip ? new IntersectionResult(clones2, clones1, clones12,
                freq2, freq1, freq21, freq12,
                clonotypes2, clonotypes1) :
                new IntersectionResult(clones1, clones2, clones12,
                        freq1, freq2, freq12, freq21,
                        clonotypes1, clonotypes2)
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
                case IntersectionType.NucleotideV:
                    if (clonotype.cdr3nt != clonotypeWrapper.clonotype.cdr3nt ||
                            clonotype.v != clonotypeWrapper.clonotype.v)
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
            if (intersectionType == IntersectionType.NucleotideV)
                return clonotype.cdr3nt.hashCode() + 31 * clonotype.v.hashCode()
            else
                return intersectionType == IntersectionType.Nucleotide ?
                        clonotype.cdr3nt.hashCode() : clonotype.cdr3aa.hashCode()
        }
    }
}