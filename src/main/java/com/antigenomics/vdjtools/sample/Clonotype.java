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

package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.ClonotypeContainer;
import com.antigenomics.vdjtools.Countable;
import com.antigenomics.vdjtools.Mutation;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * A class holding comprehensive info on a T- or B-cell clonotype.
 * CDR stands for Complementarity Determining Region
 */
public class Clonotype implements Comparable<Clonotype>, Countable {
    private ClonotypeContainer parent;
    private int count;
    private double freq;
    private final String key;

    private final int[] segmPoints;
    private final String v, d, j;
    private final String cdr3nt, cdr3aa;

    private final boolean inFrame, isComplete, noStop;


    /**
     * Creates a new clonotype explicitly setting all fields
     *
     * @param parent     parent sample
     * @param count      clonotype count, number of reads associated with this clonotype
     * @param freq       clonotype frequency, fraction of reads associated with this clonotype in the parent sample
     * @param segmPoints an array containing Variable segment end, Diversity segment start, Diversity segment end and Joining segment start in CDR3 coordinates
     * @param v          Variable segment identifier
     * @param d          Diversity segment identifier
     * @param j          Joining segment identifier
     * @param cdr3nt     nucleotide sequence of CDR3
     * @param cdr3aa     amino acid sequence of CDR3
     * @param inFrame    tells if this clonotype is in frame, i.e. it's sequence doesn't contain frameshifts
     * @param noStop     tells if the sequence of this clonotype doesn't contain any stop codon
     * @param isComplete tells if this clonotype is complete, i.e. CDR3 sequence for this clonotype is completely determined
     */
    public Clonotype(ClonotypeContainer parent, int count, double freq,
                     int[] segmPoints, String v, String d, String j,
                     String cdr3nt, String cdr3aa,
                     boolean inFrame, boolean noStop, boolean isComplete) {
        this.parent = parent;
        this.count = count;
        this.freq = freq;
        this.segmPoints = segmPoints;
        this.v = v;
        this.d = d;
        this.j = j;
        this.cdr3nt = cdr3nt;
        this.cdr3aa = cdr3aa;
        this.inFrame = inFrame;
        this.isComplete = isComplete;
        this.noStop = noStop;

        StringBuilder key = new StringBuilder(v).append(KEY_SEP).
                append(cdr3nt).append(KEY_SEP).
                append(j).append(KEY_SEP);

        key.setLength(key.length() - 1);

        this.key = key.toString();
    }

    /**
     * A copying constructor
     *
     * @param toCopy clonotype to copy all fields from
     */
    public Clonotype(Clonotype toCopy) {
        this(toCopy, toCopy.parent, toCopy.count);
    }

    /**
     * A copying constructor
     *
     * @param toCopy    clonotype to copy all fields from
     * @param newParent sample this clonotype will be assigned to
     */
    public Clonotype(Clonotype toCopy, ClonotypeContainer newParent) {
        this(toCopy, newParent, toCopy.count);
    }

    /**
     * A copying constructor
     *
     * @param toCopy    clonotype to copy all fields from
     * @param newParent sample this clonotype will be assigned to
     * @param newCount  number of reads that will be associated with this clonotype in new sample
     */
    public Clonotype(Clonotype toCopy, ClonotypeContainer newParent, int newCount) {
        this(newParent, newCount, toCopy.freq,
                toCopy.segmPoints, toCopy.v, toCopy.d, toCopy.j,
                toCopy.cdr3nt, toCopy.cdr3aa,
                toCopy.inFrame, toCopy.noStop, toCopy.isComplete);
    }

    /**
     * Gets the clonotype count
     *
     * @return number of reads that are associated with this clonotype
     */
    public int getCount() {
        return count;
    }

    /**
     * Gets the frequency of this clonotype as specified in plain-text input file
     *
     * @return fraction of reads that are associated with this clonotype according to plain-text input file
     */
    public double getFreqAsInInput() {
        return freq;
    }

    /**
     * Gets the frequency of this clonotype in its parent sample
     *
     * @return the fraction of reads that are associated with this clonotype in its parent sample
     */
    public double getFreq() {
        return count / (double) parent.getCount();
    }

    /**
     * Gets the Variable segment of this clonotype
     *
     * @return string identifier of Variable segment
     */
    public String getV() {
        return v;
    }

    /**
     * Gets the Diversity segment of this clonotype
     *
     * @return string identifier of Diversity segment
     */
    public String getD() {
        return d;
    }

    /**
     * Gets the Joining segment of this clonotype
     *
     * @return string identifier of Joining segment
     */
    public String getJ() {
        return j;
    }

    /**
     * Gets the nucleotide sequence of CDR3 region of this clonotype
     *
     * @return
     */
    public String getCdr3nt() {
        return cdr3nt;
    }

    /**
     * Gets the amino acid sequence of CDR3 region of this clonotype
     *
     * @return
     */
    public String getCdr3aa() {
        return cdr3aa;
    }

    /**
     * Returns {@code true} if this clonotype is in frame, {@code false} otherwise.
     * In case Variable segment sequence is not fully known, this usually involves checking if
     * Variable and Joining segment parts of CDR3 sequence are in-frame
     *
     * @return {@code true} if the clonotype sequence doesn't contain a frameshift, {@code false} otherwise
     */
    public boolean isInFrame() {
        return inFrame;
    }

    /**
     * Returns {@code true} if this clonotype is complete, {@code false} otherwise
     *
     * @return {@code true} if whole CDR3 sequence is known, {@code false} otherwise
     */
    public boolean isComplete() {
        return isComplete;
    }

    /**
     * Returns {@code true} if this clonotype sequence doesn't contain a stop codon, {@code false} otherwise
     *
     * @return {@code true} if the clonotype sequence doesn't contain a stop codon, {@code false} otherwise
     */
    public boolean isNoStop() {
        return noStop;
    }

    /**
     * Returns {@code true} if this clonotype sequence is coding, {@code false} otherwise
     *
     * @return {@code true} if the clonotype sequence doesn't contain a stop codon or a frameshift, {@code false} otherwise
     */
    public boolean isCoding() {
        return noStop && inFrame;
    }

    /**
     * Gets the position of the last nucleotide of Variable segment in CDR3 coordinates
     *
     * @return zero-based position of last Variable segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getVEnd() {
        return segmPoints[0];
    }

    /**
     * Gets the position of the first nucleotide of Diversity segment in CDR3 coordinates
     *
     * @return zero-based position of first Diversity segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getDStart() {
        return segmPoints[1];
    }

    /**
     * Gets the position of the last nucleotide of Diversity segment in CDR3 coordinates
     *
     * @return zero-based position of last Diversity segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getDEnd() {
        return segmPoints[2];
    }

    /**
     * Gets the position of the first nucleotide of Joining segment in CDR3 coordinates
     *
     * @return zero-based position of first Joining segment nucleotide, relative to first CDR3 nucleotide (first bp of conserved Cys codon)
     */
    public int getJStart() {
        return segmPoints[3];
    }

    /**
     * Gets the length of CDR3 region.
     *
     * @return length of CDR3 nucleotide sequence.
     */
    public int getCdr3Length() {
        return cdr3nt.length();
    }

    /**
     * Gets the number of nucleotides added between the Variable segment end and Diversity segment start
     *
     * @return number of nucleotides inserted into Variable-Diversity segment junction, or {@code -1} if Diversity segment is not defined
     */
    public int getVDIns() {
        return (segmPoints[0] >= 0 && segmPoints[1] >= 0) ? segmPoints[1] - segmPoints[0] - 1 : -1;
    }

    /**
     * Gets the number of nucleotides added between the Diversity segment end and Joining segment start
     *
     * @return number of nucleotides inserted into Diversity-Joining segment junction, or {@code -1} if Diversity segment is not defined
     */
    public int getDJIns() {
        return (segmPoints[3] >= 0 && segmPoints[2] >= 0) ? segmPoints[3] - segmPoints[2] - 1 : -1;
    }

    /**
     * Gets the total number of added nucleotides
     *
     * @return number of nucleotides inserted into Variable-Diversity and Diversity-Joining segment junctions, or {@code -1} if Diversity segment is not defined
     */
    public int getInsertSize() {
        return (segmPoints[0] >= 0 && segmPoints[1] >= 0 && segmPoints[2] >= 0 && segmPoints[3] >= 0) ? getVDIns() + getDJIns() : -1;
    }

    /**
     * Gets the number of nucleotides that are between Variable and Joining segments
     *
     * @return number of nucleotides between Variable and Joining segments, including those of Diversity segment (if present)
     */
    public int getNDNSize() {
        return (segmPoints[0] >= 0 && segmPoints[3] >= 0) ? segmPoints[3] - segmPoints[0] - 1 : 0;
    }

    /**
     * INTERNAL used for output
     *
     * @return
     */
    public String getBlank() {
        return ".";
    }

    /**
     * Gets the parent sample
     *
     * @return parent clonotype container
     */
    public ClonotypeContainer getParent() {
        return parent;
    }

    void append(Clonotype other) {
        this.count += other.count;
        this.freq += other.freq;
    }

    void recalculateFrequency() {
        this.freq = getFreq();
    }

    private static final String KEY_SEP = "_",
            MUT_SEP = "|";

    /**
     * Get unique key that is associated with a given clonotype
     * todo: consider replacing / on-demand computation
     *
     * @return clonotype key generated based on Variable segment, Joining segment, CDR3 nucleotide sequence and hypermutations
     */
    public String getKey() {
        return key;
    }

    @Override
    public int compareTo(Clonotype o) {
        return -Integer.compare(this.count, o.count);
    }

    @Override
    public boolean equals(Object o) {
        Clonotype that = (Clonotype) o;

        return key.equals(that.key) && parent.equals(that.parent);
    }

    @Override
    public int hashCode() {
        return 31 * parent.hashCode() + key.hashCode();
    }
}
