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
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.util.CommonUtil
import groovyx.gpars.GParsPool
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import sun.reflect.generics.reflectiveObjects.NotImplementedException

import java.util.concurrent.atomic.AtomicInteger

/**
 * A helper class for performing all sample intersection procedures in VDJtools
 */
class IntersectionUtil {
    private final int nPerms = 1000
    private final Random rnd = new Random(2106803L)
    private final IntersectionType intersectionType
    // todo: verbosity

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
     * Performs a pairwise intersection of sample collection
     * @param sampleCollection sample intersection to
     * @param storeIntersectedList if true will store list of intersected clonotypes from both samples (mem-consuming)
     * @param computeComplexMeasures if true will compute complex measures such as correlation and overlap frequency P-values
     * @return a list of paired intersections
     */
    PairedIntersectionMatrix intersectWithinCollection(SampleCollection sampleCollection,
                                                       boolean storeIntersectedList,
                                                       boolean computeComplexMeasures) {
        def pairs = sampleCollection.listPairs()
        def counter = new AtomicInteger()

        Collection<PairedIntersection> results = null

        GParsPool.withPool CommonUtil.THREADS, {
            results = (Collection<PairedIntersection>) pairs.collectParallel { SamplePair pair ->
                def pairedIntersection = generatePairedIntersection(pair, storeIntersectedList, computeComplexMeasures)

                println "[${new Date()} BatchIntersection] " +
                        "Intersected ${counter.incrementAndGet()} of ${pairs.size()} so far\n" +
                        "Last result\n${PairedIntersection.HEADER}\n$pairedIntersection"

                pairedIntersection
            }
        }

        new PairedIntersectionMatrix(sampleCollection, results.sort { it.parent.j }.sort { it.parent.i }, this)
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
     * Internal - collapse sample according to intersection type
     */
    private Map<ClonotypeHashWrapper, ClonotypeHashWrapper> createWrappedSample(Sample sample) {
        def wrappedSample = new HashMap<ClonotypeHashWrapper, ClonotypeHashWrapper>()

        sample.each { Clonotype clonotype ->
            def cw = new ClonotypeHashWrapper(clonotype)
            def ex = wrappedSample[cw]
            if (ex) {
                ex.freq += cw.freq
                ex.count += cw.count
                ex.clonotypes.add(cw.clonotype)
            } else {
                wrappedSample.put(cw, cw)
            }
        }

        wrappedSample
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
        def sample1 = samplePair[0], sample2 = samplePair[1]

        // Here we'll operate wrapped clonotypes and perform intersection via hash map
        // Collapse samples to hash sets

        def wrappedSample1 = createWrappedSample(sample1),
            wrappedSample2 = createWrappedSample(sample2)

        // Flip samples for speedup

        boolean flip = wrappedSample1.size() < wrappedSample2.size()

        // Overlap samples using hash-search
        // also accumulate some statistics

        (wrappedSample1, wrappedSample2) = flip ? [wrappedSample2, wrappedSample1] : [wrappedSample1, wrappedSample2]

        double freq12 = 0, freq21 = 0
        int count12 = 0, count21 = 0, clones12 = 0

        List<Clonotype> clonotypes12 = new ArrayList<>(), clonotypes21 = new ArrayList<>()

        final List<Double> x = new ArrayList<>(), y = new ArrayList<>()

        // Use larger sample as hash map and iterate through the smaller sample

        wrappedSample2.keySet().each { ClonotypeHashWrapper it ->
            def other = wrappedSample1[it]

            if (other != null) {
                clones12++

                count21 += it.count
                freq21 += it.freq

                count12 += other.count
                freq12 += other.freq

                if (computeComplexMeasures) {
                    //x.add(Math.log10(it.freq))
                    //y.add(Math.log10(other.freq))
                    x.add(it.freq)
                    y.add(other.freq)
                }

                if (storeIntersectedList) {
                    clonotypes12.addAll(other.clonotypes)
                    clonotypes21.add(it.clonotype)
                }
            }
        }

        // Flip back

        if (flip)
            (count12, count21, freq12, freq21, clonotypes12, clonotypes21) =
                    [count21, count12, freq21, freq12, clonotypes21, clonotypes12]

        // Pre-compute some complex measures

        double r, vJSD,
               freq12e, freq21e,
               freq12p, freq21p

        if (computeComplexMeasures) {
            // Correlation within intersecting set
            r = x.size() > 2 ? new PearsonsCorrelation().correlation(x as double[], y as double[]) : Double.NaN

            // Compute expected values and P-value for f12
            final List<Clonotype> clonotypes1 = new ArrayList<>(sample1.clonotypes),
                                  clonotypes2 = new ArrayList<>(sample2.clonotypes)

            // Jensen-Shannon distance
            Map<String, Double> vUsage1 = new HashMap(), vUsage2 = new HashMap<>()
            clonotypes1.each {
                vUsage1.put(it.v, (vUsage1[it.v] ?: 0d) + it.freq)
            }
            clonotypes2.each {
                vUsage2.put(it.v, (vUsage2[it.v] ?: 0d) + it.freq)
            }

            def vUsageSum1 = (double) vUsage1.values().collect().sum(),
                vUsageSum2 = (double) vUsage2.values().collect().sum()

            vJSD = [vUsage1.keySet(), vUsage2.keySet()].flatten().collect {
                double p = (vUsage1[it] ?: 0d) / vUsageSum1, q = (vUsage2[it] ?: 0d) / vUsageSum2,
                       m = (p + q) / 2.0

                (p > 0 ? (Math.log(p / m) * p) : 0d) + (q > 0 ? (Math.log(q / m) * q) : 0d)
            }.sum() / 2.0 / Math.log(2.0)

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
            vJSD = Double.NaN
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
                r, vJSD,
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
        final Clonotype coreClonotype
        public final List<Clonotype> clonotypes = new LinkedList<>()
        public double freq
        public int count

        ClonotypeHashWrapper(Clonotype coreClonotype) {
            this.coreClonotype = coreClonotype
            clonotypes.add(coreClonotype)
            this.freq = coreClonotype.freq
            this.count = coreClonotype.count
        }

        public Clonotype getClonotype() {
            coreClonotype
        }

        boolean equals(o) {
            if (this.is(o)) return true

            ClonotypeHashWrapper clonotypeWrapper = (ClonotypeHashWrapper) o

            switch (intersectionType) {
                case IntersectionType.Nucleotide:
                    if (coreClonotype.cdr3nt != clonotypeWrapper.coreClonotype.cdr3nt)
                        return false
                    break
                case IntersectionType.NucleotideV:
                    if (coreClonotype.cdr3nt != clonotypeWrapper.coreClonotype.cdr3nt ||
                            coreClonotype.v != clonotypeWrapper.coreClonotype.v)
                        return false
                    break
                case IntersectionType.AminoAcid:
                    if (coreClonotype.cdr3aa != clonotypeWrapper.coreClonotype.cdr3aa)
                        return false
                    break
                case IntersectionType.AminoAcidNonNucleotide:
                    if (coreClonotype.cdr3aa != clonotypeWrapper.coreClonotype.cdr3aa ||
                            coreClonotype.cdr3nt == clonotypeWrapper.coreClonotype.cdr3nt)
                        return false
                    break
            }

            return true
        }

        int hashCode() {
            if (intersectionType == IntersectionType.NucleotideV)
                return coreClonotype.cdr3nt.hashCode() + 31 * coreClonotype.v.hashCode()
            else
                return intersectionType == IntersectionType.Nucleotide ?
                        coreClonotype.cdr3nt.hashCode() : coreClonotype.cdr3aa.hashCode()
        }
    }
}
