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

import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample

/**
 * Class that represents spectratype, i.e. distribution of clonotypes / clonotype frequency by 
 * CDR3 amino acid / nucleotide sequence length
 */
public class Spectratype extends Histogram {
    protected final boolean aminoAcid

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
        super(unweighted, aminoAcid ? 7 : 21, aminoAcid ? 23 : 69)
        this.aminoAcid = aminoAcid
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

    @Override
    protected int getValue(Clonotype clonotype) {
        aminoAcid ? clonotype.cdr3nt.length() / 3 : clonotype.cdr3nt.length()
    }
}
