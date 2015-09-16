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

import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.sample.Sample;

import java.util.*;

public class JointSample implements Iterable<JointClonotype> {
    private final Sample[] samples;
    private final double[] intersectionFreq;
    private final double[][] intersectionFreqMatrix;
    private final long[] intersectionCount;
    private final long[][] intersectionCountMatrix;
    private final int[] intersectionDiv;
    private final int[][] intersectionDivMatrix;
    private final List<JointClonotype> jointClonotypes;
    private final double totalMeanFreq, minMeanFreq;
    private final int numberOfSamples;
    private final int count;
    private final OverlapType overlapType;
    private final boolean reverse;

    private JointSample(Sample[] samples,
                        double[] intersectionFreq, double[][] intersectionFreqMatrix,
                        long[] intersectionCount, long[][] intersectionCountMatrix,
                        int[] intersectionDiv, int[][] intersectionDivMatrix,
                        List<JointClonotype> jointClonotypes,
                        double totalMeanFreq, double minMeanFreq,
                        int numberOfSamples, int count,
                        OverlapType overlapType, boolean reverse) {
        this.samples = samples;
        this.intersectionFreq = intersectionFreq;
        this.intersectionFreqMatrix = intersectionFreqMatrix;
        this.intersectionCount = intersectionCount;
        this.intersectionCountMatrix = intersectionCountMatrix;
        this.intersectionDiv = intersectionDiv;
        this.intersectionDivMatrix = intersectionDivMatrix;
        this.jointClonotypes = jointClonotypes;
        this.totalMeanFreq = totalMeanFreq;
        this.minMeanFreq = minMeanFreq;
        this.numberOfSamples = numberOfSamples;
        this.count = count;
        this.overlapType = overlapType;
        this.reverse = reverse;
    }

    public JointSample(OverlapType overlapType, Sample[] samples) {
        this(overlapType, samples, new OccurenceJoinFilter());
    }

    public JointSample(OverlapType overlapType, Sample[] samples, JoinFilter joinFilter) {
        this.numberOfSamples = samples.length;
        this.samples = samples;
        this.intersectionDiv = new int[numberOfSamples];
        this.intersectionFreq = new double[numberOfSamples];
        this.intersectionFreqMatrix = new double[numberOfSamples][numberOfSamples];
        this.intersectionCount = new long[numberOfSamples];
        this.intersectionCountMatrix = new long[numberOfSamples][numberOfSamples];
        this.intersectionDivMatrix = new int[numberOfSamples][numberOfSamples];
        this.overlapType = overlapType;
        this.reverse = false;

        ClonotypeKeyGen clonotypeKeyGen = new ClonotypeKeyGen(overlapType);

        Map<ClonotypeKey, JointClonotype> clonotypeMap = new HashMap<>();
        int sampleIndex = 0;
        for (Sample sample : samples) {
            for (Clonotype clonotype : sample) {
                ClonotypeKey key = clonotypeKeyGen.generateKey(clonotype);

                JointClonotype jointClonotype = clonotypeMap.get(key);

                if (jointClonotype == null) {
                    clonotypeMap.put(key, jointClonotype = new JointClonotype(this));
                }

                jointClonotype.addVariant(clonotype, sampleIndex);
            }
            sampleIndex++;
        }

        this.jointClonotypes = new ArrayList<>(clonotypeMap.size() / 2);

        double totalMeanFreq = 0, minMeanFreq = 1;
        int count = 0;
        for (JointClonotype jointClonotype : clonotypeMap.values()) {
            if (joinFilter.pass(jointClonotype)) {
                jointClonotypes.add(jointClonotype);

                double meanFreq = jointClonotype.getBaseFreq();
                totalMeanFreq += meanFreq;
                minMeanFreq = Math.min(minMeanFreq, meanFreq);

                for (int i = 0; i < numberOfSamples; i++) {
                    if (jointClonotype.present(i)) {
                        count += jointClonotype.getCount(i);

                        double freq1 = jointClonotype.getFreq(i);
                        int count1 = jointClonotype.getCount(i);

                        intersectionCount[i] += count1;
                        intersectionFreq[i] += freq1;
                        intersectionDiv[i]++;

                        for (int j = i + 1; j < numberOfSamples; j++) {
                            if (jointClonotype.present(j)) {
                                double freq2 = jointClonotype.getFreq(j);
                                int count2 = jointClonotype.getCount(j);

                                intersectionFreqMatrix[i][j] += freq1;
                                intersectionFreqMatrix[j][i] += freq2;
                                intersectionCountMatrix[i][j] += count1;
                                intersectionCountMatrix[j][i] += count2;
                                intersectionDivMatrix[i][j]++;
                            }
                        }
                    }
                }
            }
        }
        this.totalMeanFreq = totalMeanFreq;
        this.minMeanFreq = minMeanFreq;
        this.count = count;

        Collections.sort(jointClonotypes);
    }

    public int getNumberOfSamples() {
        return numberOfSamples;
    }

    public Sample getSample(int sampleIndex) {
        return samples[getIndex(sampleIndex)];
    }

    public double getFreq() {
        return 1.0;
    }

    public int getDiversity() {
        return jointClonotypes.size();
    }

    public int getCount() {
        return count;
    }

    public JointClonotype getAt(int index) {
        JointClonotype jointClonotype = jointClonotypes.get(index);
        return reverse ? jointClonotype.changeParent(this) : jointClonotype;
    }

    public int getIntersectionDiv(int sampleIndex) {
        return intersectionDiv[getIndex(sampleIndex)];
    }

    public int getIntersectionDiv(int sampleIndex1, int sampleIndex2) {
        return sampleIndex1 < sampleIndex2 ? intersectionDivMatrix[sampleIndex1][sampleIndex2] :
                intersectionDivMatrix[sampleIndex2][sampleIndex1];
    }

    public double getIntersectionFreq(int sampleIndex) {
        return intersectionFreq[getIndex(sampleIndex)];
    }

    public double getIntersectionFreq(int sampleIndex1, int sampleIndex2) {
        return intersectionFreqMatrix[getIndex(sampleIndex1)][getIndex(sampleIndex2)];
    }

    public long getIntersectionCount(int sampleIndex) {
        return intersectionCount[getIndex(sampleIndex)];
    }

    public long getIntersectionCount(int sampleIndex1, int sampleIndex2) {
        return intersectionCountMatrix[getIndex(sampleIndex1)][getIndex(sampleIndex2)];
    }

    public OverlapType getOverlapType() {
        return overlapType;
    }

    public double getTotalMeanFreq() {
        return totalMeanFreq;
    }

    public double calcFreq(double meanFreq) {
        return meanFreq / totalMeanFreq;
    }

    public int calcCount(double freq) {
        return (int) (freq / minMeanFreq);
    }

    int getIndex(int sampleIndex) {
        return reverse ? (numberOfSamples - sampleIndex - 1) : sampleIndex;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        JointSample that = (JointSample) o;

        return Arrays.equals(samples, that.samples) && reverse == that.reverse;
    }

    @Override
    public int hashCode() {
        return 31 * Arrays.hashCode(samples) + (reverse ? 1231 : 1237);
    }

    @Override
    public Iterator<JointClonotype> iterator() {
        final Iterator<JointClonotype> iterator = jointClonotypes.iterator();

        if (reverse) {
            return new Iterator<JointClonotype>() {
                @Override
                public boolean hasNext() {
                    return iterator.hasNext();
                }

                @Override
                public JointClonotype next() {
                    return iterator.next().changeParent(JointSample.this);
                }

                @Override
                public void remove() {
                    throw new UnsupportedOperationException("remove");
                }
            };
        }

        return iterator;
    }

    public JointSample getReverse() throws Exception {
        if (reverse)
            throw new Exception("Already reversed, for the sake of performance multiple reverse is disabled");

        return new JointSample(samples,
                intersectionFreq, intersectionFreqMatrix,
                intersectionCount, intersectionCountMatrix,
                intersectionDiv, intersectionDivMatrix,
                jointClonotypes,
                totalMeanFreq, minMeanFreq,
                numberOfSamples, count, overlapType, true);
    }
}
