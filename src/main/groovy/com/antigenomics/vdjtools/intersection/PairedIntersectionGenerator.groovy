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
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SamplePair

class PairedIntersectionGenerator {
    private final Collection<ClonotypeWrapper> wrappedSample1, wrappedSample2
    private final Sample sample1, sample2

    final IntersectionType intersectionType

    PairedIntersectionGenerator(SamplePair samplePair,
                                IntersectionType intersectionType) {
        this(samplePair.sample1, samplePair.sample2, intersectionType)
    }

    PairedIntersectionGenerator(Sample sample1, Sample sample2,
                                IntersectionType intersectionType) {
        this.intersectionType = intersectionType
        this.sample1 = sample1
        this.sample2 = sample2
        this.wrappedSample1 = sample1.clonotypes.collect { new ClonotypeWrapper(it, intersectionType) }
        this.wrappedSample2 = sample2.clonotypes.collect { new ClonotypeWrapper(it, intersectionType) }
    }

    PairedIntersectionResult intersect(boolean store) {
        def intersection = new HashMap<ClonotypeWrapper, ClonotypeWrapper>()

        // flip is just for speedup
        boolean flip = wrappedSample1.size() < wrappedSample2.size()

        Collection<ClonotypeWrapper> _sample1, _sample2

        (_sample1, _sample2) = flip ? [wrappedSample2, wrappedSample1] : [wrappedSample1, wrappedSample2]

        double freq12 = 0, freq21 = 0
        int count12 = 0, count21 = 0, clones12 = 0

        _sample1.each {
            intersection.put(it, it)
        }

        final List<Clonotype> clonotypes12 = new ArrayList<>(),
                              clonotypes21 = new ArrayList<>()

        _sample2.each {
            def other = intersection[it]

            if (other != null) {
                clones12++

                count12 += it.clonotype.count
                freq12 += it.clonotype.freq

                count21 += other.clonotype.count
                freq21 += other.clonotype.freq

                if (store) {
                    clonotypes12.add(other.clonotype)
                    clonotypes21.add(it.clonotype)
                }
            }
        }

        flip ? new PairedIntersectionResult(
                sample1, sample2,
                clones12,
                count21, count12,
                freq21, freq12,
                clonotypes21, clonotypes12)
                :
                new PairedIntersectionResult(
                        sample1, sample2,
                        clones12,
                        count12, count21,
                        freq12, freq21,
                        clonotypes12, clonotypes21)
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