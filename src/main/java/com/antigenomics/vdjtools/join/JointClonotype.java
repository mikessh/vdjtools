/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.ClonotypeWrapper;
import com.antigenomics.vdjtools.Misc;
import com.antigenomics.vdjtools.sample.Clonotype;
import org.apache.commons.math3.stat.inference.ChiSquareTest;

import java.util.LinkedList;
import java.util.List;

public class JointClonotype implements Comparable<JointClonotype>, ClonotypeWrapper {
    private final JointSample parent;
    private final List[] variantsBySample;
    private final int[] counts;
    private int peak = -1, occurences = -1;
    private Clonotype representative = null;
    private double meanFreq = -1;
    protected double samplingPValue = 1;

    private JointClonotype(JointSample parent,
                           List[] variantsBySample, int[] counts,
                           int peak, Clonotype representative, double meanFreq) {
        this.parent = parent;
        this.variantsBySample = variantsBySample;
        this.counts = counts;
        this.peak = peak;
        this.representative = representative;
        this.meanFreq = meanFreq;
    }

    public JointClonotype(JointSample parent) {
        this.parent = parent;
        this.variantsBySample = new List[parent.getNumberOfSamples()];
        this.counts = new int[parent.getNumberOfSamples()];
    }

    void addVariant(Clonotype variant, int sampleIndex) {
        List<Clonotype> variants = variantsBySample[sampleIndex];
        if (variants == null)
            variantsBySample[sampleIndex] = (variants = new LinkedList());
        variants.add(variant);
        counts[sampleIndex] += variant.getCount();
    }

    public JointSample getParent() {
        return parent;
    }

    public int getPeak() {
        if (peak < 0) {
            int peakCount = -1;
            for (int i = 0; i < counts.length; i++) {
                if (counts[i] > peakCount) {
                    peak = i;
                    peakCount = counts[i];
                }
            }
        }
        return parent.getIndex(peak);
    }

    public int getOccurences() {
        if (occurences < 0) {
            occurences = 0;
            for (int i = 0; i < counts.length; i++) {
                if (present(i))
                    occurences++;
            }
        }
        return occurences;
    }

    public double getSamplingPValue() {
        return samplingPValue;
    }

    @Override
    public Clonotype getClonotype() {
        if (representative == null) {
            int max = 0;
            getPeak();
            for (Object o : variantsBySample[peak]) {
                Clonotype clonotype = (Clonotype) o;
                if (clonotype.getCount() > max) {
                    max = clonotype.getCount();
                    representative = clonotype;
                }
            }
        }
        return representative;
    }

    public int getNumberOfVariants(int sampleIndex) {
        return variantsBySample[parent.getIndex(sampleIndex)].size();
    }

    public int getCount(int sampleIndex) {
        return counts[parent.getIndex(sampleIndex)];
    }

    /**
     * Gets clonotype frequency in a given sample
     *
     * @param sampleIndex
     * @return
     */
    public double getFreq(int sampleIndex) {
        return getCount(sampleIndex) / (double) parent.getCount(sampleIndex);
    }

    /**
     * Gets clonotype frequency at a given sample relative to overlap size of this sample
     *
     * @param sampleIndex
     * @return
     */
    public double getFreqWithinIntersection(int sampleIndex) {
        return getFreq(sampleIndex) / parent.getIntersectionFreq(sampleIndex);
    }

    /**
     * Gets a representative (geometric mean) frequency of a given joint clonotype
     *
     * @return
     */
    public double getBaseFreq() {
        if (meanFreq < 0) {
            meanFreq = 1;
            for (int i = 0; i < parent.getNumberOfSamples(); i++) {
                meanFreq *= (getFreq(i) + Misc.JITTER);
            }
            meanFreq = Math.pow(meanFreq, 1.0 / (double) parent.getNumberOfSamples());
        }
        return meanFreq;
    }

    /**
     * Gets a representative (geometric mean) frequency of a given joint clonotype.
     * The frequency is normalized so that the total mean frequency will sum to 1.0 in JointSample,
     * thus allowing to use it during sample output
     *
     * @return
     */
    public double getFreq() {
        return parent.calcFreq(getBaseFreq());
    }

    /**
     * Gets a representative count of a given joint clonotype, which is proportional to representative
     * (geometric mean) frequency. The count is scaled so it is equal to 1 for the smallest JointClonotype,
     * thus allowing to use it during sample output
     *
     * @return
     */
    public int getCount() {
        return parent.calcCount(getBaseFreq());
    }

    public boolean present(int sampleIndex) {
        return counts[parent.getIndex(sampleIndex)] > 0;
    }

    JointClonotype changeParent(JointSample newParent) {
        return new JointClonotype(newParent, variantsBySample, counts, peak, representative, meanFreq);
    }

    @Override
    public int compareTo(JointClonotype o) {
        return Double.compare(o.getBaseFreq(), this.getBaseFreq());
    }
}
