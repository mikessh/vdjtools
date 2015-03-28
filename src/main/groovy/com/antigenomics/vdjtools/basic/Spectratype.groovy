/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.sample.Sample

/**
 * Class that represents spectratype, i.e. distribution of clonotypes / clonotype frequency by 
 * CDR3 amino acid / nucleotide sequence length
 */
public class Spectratype {
    private final boolean aminoAcid, unweighted
    private final int min, max, span

    private final double[] innerHist
    private final int[] lengths

    private int count = 0
    private double freq = 0

    /**
     * Creates a spectratype instance with a given sample. All calculations are performed within constructor.
     * @param sample sample to initialize with
     * @param aminoAcid will use amino acid CDR3 sequence if true. Will use nucleotide sequence otherwise
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public Spectratype(Sample sample, boolean aminoAcid, boolean unweighted) {
        this(aminoAcid, unweighted)
        addAll(sample)
    }

    /**
     * Creates a spectratype instance with a given sample. All calculations are performed within constructor.
     * Internal constructor 
     * @param sample sample sample to initialize with
     * @param intersectionType overlap type to deduce whether amino acid or nucleotide sequence CDR3 should be used
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public Spectratype(Sample sample, OverlapType intersectionType, boolean unweighted) {
        this(intersectionType.aminoAcid, unweighted)
        addAll(sample)
    }

    /**
     * Creates a blank spectratype instance.
     * @param intersectionType overlap type to deduce whether amino acid or nucleotide sequence CDR3 should be used
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
    public Spectratype(OverlapType intersectionType, boolean unweighted) {
        this(intersectionType.aminoAcid, unweighted)
    }

    /**
     * Creates a blank spectratype instance. 
     * @param aminoAcid will use amino acid CDR3 sequence if true. Will use nucleotide sequence otherwise
     * @param unweighted will count each unique clonotype once if set to true. Will weight each clonotype by its frequency otherwise
     */
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

        this.span = max - min + 1
        this.innerHist = new double[span]
        this.lengths = (min..max) as double[]
        this.HEADER = lengths.collect().join("\t")
    }

    /**
     * INTERNAL
     */
    private int bin(Clonotype clonotype) {
        Math.min(max, Math.max(min, (aminoAcid ? getCdr3AaLen(clonotype) : clonotype.cdr3nt.length()))) - min
    }

    /**
     * INTERNAL
     */
    private static int getCdr3AaLen(Clonotype clonotype) {
        (int) (clonotype.inFrame ? clonotype.cdr3aa.length() : (clonotype.cdr3nt.length() / 3))
    }

    /**
     * Update spectratype with a single clonotype 
     * @param clonotype clonotype to add
     */
    public void add(Clonotype clonotype) {
        if (unweighted) {
            this.innerHist[bin(clonotype)]++
            freq++
        } else {
            this.innerHist[bin(clonotype)] += clonotype.freq
            freq += clonotype.freq
        }
        count++
    }

    /**
     * Update spectratype with a set of clonotypes 
     * @param sample sample to add
     */
    public void addAll(Iterable<Clonotype> sample) {
        sample.each { Clonotype clonotype ->
            add(clonotype)
        }
    }

    /**
     * Combines two spectratypes 
     * @param other spectratype to update with
     */
    public void addAll(Spectratype other) {
        if (other.aminoAcid != this.aminoAcid || other.unweighted != this.unweighted)
            throw new Exception("Can't add Spectratype of different type")

        this.count += other.count
        this.freq += other.freq
        other.innerHist.eachWithIndex { it, ind ->
            innerHist[ind] += it
        }
    }

    /**
     * Update spectratype with a set of clonotypes, but omits several top clonotypes.
     * Internal, used for fancy (detalized) spectratype.
     * @param sample sample to add
     * @param top top clonotypes which will be omitted, but rather returned as list
     * @return
     */
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

    /**
     * Clears all spectratype counters 
     */
    public void clear() {
        for (int i = 0; i < span; i++)
            this.innerHist[i] = 0
        count = 0
        freq = 0
    }

    /**
     * Gets the spectratype histogram. Values correspond to clonotype CDR3 lengths from {@code min} to {@code max}
     * @return
     */
    public double[] getHistogram() {
        getHistogram(true)
    }

    /**
     * Gets the spectratype histogram. Values correspond to clonotype CDR3 lengths from {@code min} to {@code max}
     * @param normalized if true will normalize clonotype to total sum of 1.
     * @return
     */
    public double[] getHistogram(boolean normalized) {
        def spectratype = new double[span]
        
        double _freq = (normalized && freq > 0) ? freq : 1.0

        for (int i = 0; i < span; i++)
            spectratype[i] = this.innerHist[i] / _freq

        spectratype
    }

    /**
     * Gets min CDR3 length of clonotype (lower boundary on spectratype length)
     * @return
     */
    int getMin() {
        min
    }

    /**
     * Gets max CDR3 length of clonotype (upper boundary on spectratype length)
     * @return
     */
    int getMax() {
        max
    }

    /**
     * Gets the number of bins in spectratype histogram, i.e. {@code max - min + 1}
     * @return
     */
    public int getSpan() {
        span
    }

    /**
     * Gets total frequency of summarized clonotypes 
     * @return
     */
    public double getFreq() {
        freq
    }

    /**
     * Gets total read count of summarized clonotypes
     * @return
     */
    public int getCount() {
        count
    }

    /**
     * An array of CDR3 lengths that correspond to spectratype bins
     */
    public int[] getLengths() {
        return lengths
    }

    /**
     * Header string, used for tabular output
     */
    public final String HEADER

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        histogram.collect().join("\t")
    }
}
