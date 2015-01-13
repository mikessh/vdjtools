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
 *
 * Last modified on 13.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.sample.Sample
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

/**
 * This class provides some basic stats that could be computed for a RepSeq sample:
 * read count, observed diversity, characteristic clone size, cdr3/N-N insert/N-D-N sequence size,
 * fraction of non-coding clonotypes and repertoire convergence 
 */
public class BasicStats {
    private final DescriptiveStatistics cloneSize, cloneSizeGeom, cdr3ntLength, insertSize, ndnSize
    private int ncDiversity
    private double ncFrequency
    private final long count
    private final int diversity
    private final double convergence

    /**
     * Create an instance of BasicStats class. All computations are done within constructor.
     * @param sample sample to be analyzed
     */
    public BasicStats(Sample sample) {
        this.count = sample.count
        this.diversity = sample.diversity

        this.cloneSize = new DescriptiveStatistics()
        this.cloneSizeGeom = new DescriptiveStatistics()
        this.cdr3ntLength = new DescriptiveStatistics()
        this.insertSize = new DescriptiveStatistics()
        this.ndnSize = new DescriptiveStatistics()

        def clonotypeKeyGen = new ClonotypeKeyGen(IntersectionType.AminoAcidVJ)
        def aaSet = new HashSet<ClonotypeKey>()

        sample.each {
            cloneSize.addValue(it.freq)
            cloneSizeGeom.addValue(Math.log10(it.freq))
            cdr3ntLength.addValue(it.cdr3nt.length())

            def x = it.insertSize
            if (x > -1)
                insertSize.addValue(x)

            ndnSize.addValue(it.NDNSize)
            if (!it.coding) {
                ncDiversity++
                ncFrequency += it.freq
            }

            aaSet.add(clonotypeKeyGen.generateKey(it))
        }

        this.convergence = diversity / (double) aaSet.size() - 1
    }

    /**
     * Gets mean clonotype frequency 
     * @return
     */
    public double getMeanFrequency() {
        cloneSize.mean
    }

    /**
     * Gets geometric mean of clonotype frequency 
     * @return
     */
    public double getGeomMeanFrequency() {
        Math.pow(10, cloneSizeGeom.mean)
    }

    /**
     * Gets mean CDR3 nucleotide sequence length 
     * @return
     */
    public double getMeanCdr3ntLength() {
        cdr3ntLength.mean
    }

    /**
     * Gets mean NDN region size, i.e. mean number of nucleotides between V segment end and J segment start 
     * @return
     */
    public double getMeanNDNSize() {
        ndnSize.mean
    }

    /**
     * Gets mean insert size, i.e. mean number of nucleotides in V-D and D-J regions if D segment is determined or
     * V-J regions if not 
     * @return
     */
    public double getMeanInsertSize() {
        insertSize.mean
    }

    /**
     * Gets the number of reads in sample 
     * @return
     */
    public long getCount() {
        count
    }

    /**
     * Gets the total number of clonotypes in sample
     * @return
     */
    public int getDiversity() {
        diversity
    }

    /**
     * Gets the total number of non-coding (containing a stop or frameshift) clonotypes in sample 
     * @return
     */
    public int getNcDiversity() {
        ncDiversity
    }

    /**
     * Gets the cumulative frequency of non-coding (containing a stop or frameshift) clonotypes in sample 
     * @return
     */
    public double getNcFrequency() {
        ncFrequency
    }

    /**
     * Gets the convergence level of sample, i.e. number of CDR3 nucleotide sequences per CDR3 amino acid sequence 
     * @return
     */
    public double getConvergence() {
        return convergence
    }

    /**
     * Header used for tabular output 
     */
    public static final String HEADER = "count\tdiversity\t" +
            "mean_frequency\tgeomean_frequency\t" +
            "nc_diversity\tnc_frequency\t" +
            "mean_cdr3nt_length\tmean_insert_size\tmean_ndn_size\t" +
            "convergence"

    /**
     * Plain text row for tabular output 
     */
    @Override
    public String toString() {
        [count, diversity,
         meanFrequency, geomMeanFrequency,
         ncDiversity, ncFrequency,
         meanCdr3ntLength, meanInsertSize, meanNDNSize,
         convergence
        ].flatten().join("\t")
    }
}
