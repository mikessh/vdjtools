/**
 * Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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

package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.Clonotype;
import com.antigenomics.vdjtools.intersection.IntersectionUtil;
import com.antigenomics.vdjtools.sample.Sample;

import java.util.*;

public class JointSample implements Iterable<JointClonotype> {
    private final Sample[] samples;
    private final double[] intersectionFreq;
    private final double[][] intersectionFreqMatrix;
    private final int[] intersectionDiv;
    private final int[][] intersectionDivMatrix;
    private final List<JointClonotype> jointClonotypes;
    private final int sampleCount, detectionThreshold;

    public JointSample(IntersectionUtil intersectionUtil, Sample[] samples) {
        this(intersectionUtil, samples, 2);
    }

    public JointSample(IntersectionUtil intersectionUtil, Sample[] samples, int detectionThreshold) {
        this.sampleCount = samples.length;
        this.samples = samples;
        this.detectionThreshold = detectionThreshold;
        this.intersectionDiv = new int[sampleCount];
        this.intersectionFreq = new double[sampleCount];
        this.intersectionFreqMatrix = new double[sampleCount][sampleCount];
        this.intersectionDivMatrix = new int[sampleCount][sampleCount];

        Map<String, JointClonotype> clonotypeMap = new HashMap<>();
        int sampleIndex = 0;
        for (Sample sample : samples) {
            for (Clonotype clonotype : sample) {
                String key = intersectionUtil.generateKey(clonotype);

                JointClonotype jointClonotype = clonotypeMap.get(key);

                if (jointClonotype == null) {
                    clonotypeMap.put(key, jointClonotype = new JointClonotype(key, this));
                }

                jointClonotype.addVariant(clonotype, sampleIndex);
            }
            sampleIndex++;
        }

        this.jointClonotypes = new ArrayList<>(clonotypeMap.size() / 2);

        for (JointClonotype jointClonotype : clonotypeMap.values()) {
            if (jointClonotype.getNumberOfSamplesWhereDetected() >= detectionThreshold) {
                jointClonotypes.add(jointClonotype);
                for (int i = 0; i < sampleCount; i++) {
                    double freq1 = jointClonotype.getFreq(i);
                    if (freq1 > 0) {
                        intersectionFreq[i] += freq1;
                        intersectionDiv[i]++;
                        for (int j = i + 1; j < sampleCount; j++) {
                            double freq2 = jointClonotype.getFreq(j);
                            if (freq2 > 0) {
                                intersectionFreqMatrix[i][j] += Math.sqrt(freq1 * freq2);
                                intersectionDivMatrix[i][j]++;
                            }
                        }
                    }
                }
            }
        }

        Collections.sort(jointClonotypes);
    }

    public int getSampleCount() {
        return sampleCount;
    }

    public Sample[] getSamples() {
        return samples;
    }

    public int size() {
        return jointClonotypes.size();
    }

    public JointClonotype getAt(int index) {
        return jointClonotypes.get(index);
    }

    public int getDetectionThreshold() {
        return detectionThreshold;
    }

    public int getIntersectionDiv(int sampleIndex) {
        return intersectionDiv[sampleIndex];
    }

    public int getIntersectionDiv(int sampleIndex1, int sampleIndex2) {
        return sampleIndex1 < sampleIndex2 ? intersectionDivMatrix[sampleIndex1][sampleIndex2] :
                intersectionDivMatrix[sampleIndex2][sampleIndex1];
    }

    public double getIntersectionFreq(int sampleIndex) {
        return intersectionFreq[sampleIndex];
    }

    public double getIntersectionFreq(int sampleIndex1, int sampleIndex2) {
        return sampleIndex1 < sampleIndex2 ? intersectionFreqMatrix[sampleIndex1][sampleIndex2] :
                intersectionFreqMatrix[sampleIndex2][sampleIndex1];
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        JointSample that = (JointSample) o;

        return Arrays.equals(samples, that.samples);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(samples);
    }

    @Override
    public Iterator<JointClonotype> iterator() {
        return jointClonotypes.iterator();
    }
}
