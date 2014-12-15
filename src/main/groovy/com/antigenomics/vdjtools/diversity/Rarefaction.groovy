/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 * Last modified on 15.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.MathUtil


class Rarefaction {
    private final Clonotype[] flattenedClonotypes
    private final int totalReads

    public Rarefaction(Sample sample) {
        if (sample.count > Integer.MAX_VALUE)
            throw new RuntimeException("Couldn't downsample samples with > ${Integer.MAX_VALUE} cells")

        this.totalReads = (int) sample.count
        this.flattenedClonotypes = new Clonotype[totalReads]

        int counter = 0
        sample.each {
            for (int i = 0; i < it.count; i++)
                flattenedClonotypes[counter++] = it
        }
    }

    public RarefactionCurve build() {
        build(IntersectionType.Strict)
    }

    public RarefactionCurve build(IntersectionType intersectionType) {
        build(intersectionType, 3, 101)
    }

    public RarefactionCurve build(IntersectionType intersectionType, int numberOfResamples, int numberOfPoints) {
        def x = new int[numberOfPoints], y = new int[numberOfPoints][numberOfResamples]
        def clonotypeKeyGen = new ClonotypeKeyGen(intersectionType)

        int step = totalReads / (numberOfPoints - 1)

        for (int k = 0; k < numberOfResamples; k++) {
            // Count unique clonotypes using hash set
            final HashSet<ClonotypeKey> countingHash = new HashSet<>()

            MathUtil.shuffle(flattenedClonotypes) // anyways, we're ordered by size initially

            // Fill hash step by step
            for (int i = 1; i < numberOfPoints - 1; i++) {
                int xx = (i - 1) * step

                for (int j = 0; j < step; j++)
                    countingHash.add(clonotypeKeyGen.generateKey(flattenedClonotypes[xx + j]))

                x[i] = i * step
                y[i][k] = countingHash.size()
            }

            // Last step
            if (k == 0) {
                for (int j = x[numberOfPoints - 1]; j < totalReads; j++)
                    countingHash.add(clonotypeKeyGen.generateKey(flattenedClonotypes[j]))

                int totalDiv = countingHash.size()
                x[numberOfPoints - 1] = totalReads

                for (int kk = 0; kk < numberOfResamples; kk++)
                    y[numberOfPoints - 1][kk] = totalDiv
            }
        }

        new RarefactionCurve(x, y)
    }
}
