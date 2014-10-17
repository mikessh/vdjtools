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

import java.util.LinkedList;
import java.util.List;

public class JointClonotype implements Comparable<JointClonotype> {
    private final String key;
    private final JointSample parent;
    private final List[] variantsBySample;
    private final int[] counts;
    private static final double JITTER = 1e-9; // upper limit on current precision of RepSeq

    public JointClonotype(String key, JointSample parent) {
        this.key = key;
        this.parent = parent;
        this.variantsBySample = new List[parent.getSampleCount()];
        this.counts = new int[parent.getSampleCount()];
    }

    public void addVariant(Clonotype variant, int sampleIndex) {
        List<Clonotype> variants = variantsBySample[sampleIndex];
        if (variants == null)
            variantsBySample[sampleIndex] = (variants = new LinkedList());
        variants.add(variant);
        counts[sampleIndex] += variant.getCount();
    }

    public int getNumberOfSamplesWhereDetected() {
        int numberOfSamplesWhereDetected = 0;

        for (int count : counts)
            if (count > 0)
                numberOfSamplesWhereDetected++;

        return numberOfSamplesWhereDetected;
    }

    public int getNumberOfVariants(int sampleIndex) {
        return variantsBySample[sampleIndex].size();
    }

    public int getCount(int sampleIndex) {
        return counts[sampleIndex];
    }

    public double getFreq(int sampleIndex) {
        return getCount(sampleIndex) / (double) parent.getSamples()[sampleIndex].getCount();
    }

    public double getMeanFreq() {
        double meanFreq = 1;
        for (int i = 0; i < parent.getSampleCount(); i++) {
            meanFreq *= (getFreq(i) + JITTER);
        }
        return Math.pow(meanFreq, 1.0 / (double) parent.getSampleCount());
    }

    @Override
    public int compareTo(JointClonotype o) {
        return -Double.compare(this.getMeanFreq(), o.getMeanFreq());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        JointClonotype that = (JointClonotype) o;

        return key.equals(that.key) && parent.equals(that.parent);
    }

    @Override
    public int hashCode() {
        return 31 * key.hashCode() + parent.hashCode();
    }
}
