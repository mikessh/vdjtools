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
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import sun.reflect.generics.reflectiveObjects.NotImplementedException

/**
 * A helper class for performing all sample intersection procedures in VDJtools
 */
class IntersectionUtil {
    private final int nPerms = 1000
    private final Random rnd = new Random(2106803L)
    private final IntersectionType intersectionType

    /**
     * Creates a helper class to perform sample intersection
     * @param intersectionType type of clonotype intersection
     */
    IntersectionUtil(IntersectionType intersectionType) {
        this.intersectionType = intersectionType
    }

    /**
     * Generates a key (e.g. TRBVx-TGTGCTAGCTGGGCATGCTTC) according to specified intersection type
     * @param clonotype clonotype to generate a key for
     * @return string key of a clonotype
     */
    String generateKey(Clonotype clonotype) {
        switch (intersectionType) {
            case IntersectionType.NucleotideV:
                return clonotype.cdr3nt + "\t" + clonotype.v
            case IntersectionType.Nucleotide:
                return clonotype.cdr3nt
            case IntersectionType.AminoAcid:
                return clonotype.cdr3aa
            default:
                throw new NotImplementedException()
        }
    }

    /**
     * Performs a paired intersection and computes several intersection measures.
     * Measures will have index 1, 2 or 12 and 21 in the resulting object.
     * E.g. freq1 will be the total frequency in first sample,
     * while freq12 and freq21 will be the total frequnecy of overlapping clonotypes in
     * first and second samples respectively.
     *
     * To be used for intersection among individual samples not within a SampleCollection
     * Will calculate and store all possible information characterizing the intersection
     * @param sample1 first sample
     * @param sample2 second sample
     * @param storeIntersectedList if true will store list of intersected clonotypes from both samples (mem-consuming)
     * @param computeComplexMeasures if true will compute complex measures such as correlation and overlap frequency P-values
     * @return an object holding information on paired intersection
     */
    PairedIntersection generatePairedIntersection(Sample sample1, Sample sample2) {
        generatePairedIntersection(sample1, sample2, true, true)
    }

    /**
     * Performs a paired intersection and computes several intersection measures.
     * Measures will have index 1, 2 or 12 and 21 in the resulting object.
     * E.g. freq1 will be the total frequency in first sample,
     * while freq12 and freq21 will be the total frequnecy of overlapping clonotypes in
     * first and second samples respectively.
     *
     * To be used for intersection among individual samples not within a SampleCollection
     * @param sample1 first sample
     * @param sample2 second sample
     * @param storeIntersectedList if true will store list of intersected clonotypes from both samples (mem-consuming)
     * @param computeComplexMeasures if true will compute complex measures such as correlation and overlap frequency P-values
     * @return an object holding information on paired intersection
     */
    PairedIntersection generatePairedIntersection(Sample sample1, Sample sample2,
                                                  boolean storeIntersectedList,
                                                  boolean computeComplexMeasures) {
        generatePairedIntersection(new SamplePair(sample1, sample2), storeIntersectedList, computeComplexMeasures)
    }

    /**
     * Performs a paired intersection and computes several intersection measures.
     * Measures will have index 1, 2 or 12 and 21 in the resulting object.
     * E.g. freq1 will be the total frequency in first sample,
     * while freq12 and freq21 will be the total frequnecy of overlapping clonotypes in
     * first and second samples respectively.
     *
     * Will calculate and store all possible information characterizing the intersection
     * @param samplePair a pair of samples to intersect
     * @return an object holding information on paired intersection
     */
    PairedIntersection generatePairedIntersection(SamplePair samplePair) {
        generatePairedIntersection(samplePair, true, true)
    }

    /**
     * Performs a paired intersection and computes several intersection measures.
     * Measures will have index 1, 2 or 12 and 21 in the resulting object.
     * E.g. freq1 will be the total frequency in first sample,
     * while freq12 and freq21 will be the total frequnecy of overlapping clonotypes in
     * first and second samples respectively.
     * @param samplePair a pair of samples to intersect
     * @param storeIntersectedList if true will store list of intersected clonotypes from both samples (mem-consuming)
     * @param computeComplexMeasures if true will compute complex measures such as correlation and overlap frequency P-values
     * @return an object holding information on paired intersection
     */
    PairedIntersection generatePairedIntersection(SamplePair samplePair,
                                                  boolean storeIntersectedList,
                                                  boolean computeComplexMeasures) {
        def sample1 = samplePair.sample1, sample2 = samplePair.sample2

        // Here we'll operate wrapped clonotypes and perform intersection via hash map

        def wrappedSample1 = sample1.clonotypes.collect { new ClonotypeHashWrapper(it) },
            wrappedSample2 = sample2.clonotypes.collect { new ClonotypeHashWrapper(it) }

        def intersection = new HashMap<ClonotypeWrapper, ClonotypeWrapper>()

        // Flip samples for speedup

        boolean flip = wrappedSample1.size() < wrappedSample2.size()

        Collection<ClonotypeWrapper> _sample1, _sample2

        // Overlap samples using hash-search
        // also accumulate some statistics

        (_sample1, _sample2) = flip ? [wrappedSample2, wrappedSample1] : [wrappedSample1, wrappedSample2]

        double freq12 = 0, freq21 = 0
        int count12 = 0, count21 = 0, clones12 = 0

        // - put larger sample into hash

        _sample1.each {
            intersection.put(it, it)
        }

        List<Clonotype> clonotypes12 = new ArrayList<>(), clonotypes21 = new ArrayList<>()

        final List<Double> x = new ArrayList<>(), y = new ArrayList<>()

        // - iterate through the smaller sample

        _sample2.each {
            def other = intersection[it]

            if (other != null) {
                clones12++

                count21 += it.clonotype.count
                freq21 += it.clonotype.freq

                count12 += other.clonotype.count
                freq12 += other.clonotype.freq

                if (computeComplexMeasures) {
                    x.add(Math.log10(it.clonotype.freq))
                    y.add(Math.log10(other.clonotype.freq))
                }

                if (storeIntersectedList) {
                    clonotypes12.add(other.clonotype)
                    clonotypes21.add(it.clonotype)
                }
            }
        }

        // Flip back

        if (flip)
            (count12, count21, freq12, freq21, clonotypes12, clonotypes21) =
                    [count21, count12, freq21, freq12, clonotypes21, clonotypes12]

        // Pre-compute some complex measures

        double r,
               freq12e, freq21e,
               freq12p, freq21p

        if (computeComplexMeasures) {
            // Correlation within intersecting set
            r = x.size() > 2 ? new PearsonsCorrelation().correlation(x as double[], y as double[]) : Double.NaN

            // Compute expected values and P-value for f12
            final List<Clonotype> clonotypes1 = new ArrayList<>(sample1.clonotypes),
                                  clonotypes2 = new ArrayList<>(sample1.clonotypes)

            Collections.shuffle(clonotypes1)
            Collections.shuffle(clonotypes2)

            freq12e = 0
            freq21e = 0
            freq12p = 0
            freq21p = 0

            // Permutations
            for (int i = 0; i < nPerms; i++) {
                int rnd1 = rnd.nextInt(clonotypes1.size() - clones12),
                    rnd2 = rnd.nextInt(clonotypes2.size() - clones12)

                // Sample frequencies
                def freq12s = (double) clonotypes1[rnd1..(rnd1 + clones12)].sum { it.freq },
                    freq21s = (double) clonotypes2[rnd2..(rnd2 + clones12)].sum { it.freq }

                // Expected frequencies
                freq12e += freq12s
                freq21e += freq21s

                // Empirical p-values
                if (freq12s > freq12)
                    freq12p += 1
                if (freq21s > freq21)
                    freq21p += 1
            }

            freq12e /= nPerms
            freq21e /= nPerms
            freq12p /= nPerms
            freq21p /= nPerms
        } else {
            r = Double.NaN
            freq12e = Double.NaN
            freq21e = Double.NaN
            freq12p = Double.NaN
            freq21p = Double.NaN
        }

        new PairedIntersection(
                samplePair,
                clones12,
                count12, count21,
                freq12, freq21,
                freq12e, freq21e,
                freq12p, freq21p,
                r,
                clonotypes12, clonotypes21)
    }

    /**
     * Wraps the clonotype to override its equals and hashCode methods
     * according to policy specified by intersectionType
     * @param clonotype clonotype to wrap
     * @return wrapped clonotype
     */
    public ClonotypeWrapper wrap(Clonotype clonotype) {
        new ClonotypeHashWrapper(clonotype)
    }

    /**
     * Internal wrapper with overridden equals and hash code for purposes of parent class
     */
    private class ClonotypeHashWrapper implements ClonotypeWrapper {
        final Clonotype clonotype

        ClonotypeHashWrapper(Clonotype clonotype) {
            this.clonotype = clonotype
        }

        boolean equals(o) {
            if (this.is(o)) return true

            ClonotypeHashWrapper clonotypeWrapper = (ClonotypeHashWrapper) o

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
