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
import com.antigenomics.vdjtools.ClonotypeWrapper
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SamplePair

class PairedIntersectionGenerator {
    private final Collection<ClonotypeWrapper> wrappedSample1, wrappedSample2
    private final Sample sample1, sample2

    final IntersectionUtil intersectionUtil

    PairedIntersectionGenerator(SamplePair samplePair,
                                IntersectionType intersectionType) {
        this(samplePair.sample1, samplePair.sample2, intersectionType)
    }

    PairedIntersectionGenerator(Sample sample1, Sample sample2,
                                IntersectionType intersectionType) {
        this.intersectionUtil = new IntersectionUtil(intersectionType)
        this.sample1 = sample1
        this.sample2 = sample2
        this.wrappedSample1 = sample1.clonotypes.collect { intersectionUtil.wrap(it) }
        this.wrappedSample2 = sample2.clonotypes.collect { intersectionUtil.wrap(it) }
    }

    PairedIntersection intersect(boolean store) {
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

        flip ? new PairedIntersection(
                sample1, sample2,
                clones12,
                count21, count12,
                freq21, freq12,
                clonotypes21, clonotypes12)
                :
                new PairedIntersection(
                        sample1, sample2,
                        clones12,
                        count12, count21,
                        freq12, freq21,
                        clonotypes12, clonotypes21)
    }
}