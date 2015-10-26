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

package com.antigenomics.vdjtools.preprocess;

import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.Sample;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.core.tree.NeighborhoodIterator;
import com.milaboratory.core.tree.SequenceTreeMap;
import com.milaboratory.core.tree.TreeSearchParameters;

import java.util.Map;
import java.util.stream.Collectors;

public class Corrector {
    private final TreeSearchParameters treeSearchParameters;
    private final double logRatioThreshold;

    public Corrector() {
        this(2, 0.05f);
    }

    public Corrector(int maxMismatches,
                     float ratioThreshold) {
        this.treeSearchParameters = new TreeSearchParameters(maxMismatches, 0, 0);
        this.logRatioThreshold = -Math.log10(ratioThreshold);
    }

    public Sample correct(Sample sample) {
        SequenceTreeMap<NucleotideSequence, Clonotype> sequenceTreeMap = new SequenceTreeMap<>(NucleotideSequence.ALPHABET);

        for (int i = sample.getDiversity() - 1; i >= 0; i--) {
            Clonotype clonotype = sample.getAt(i);
            sequenceTreeMap.put(clonotype.getCdr3ntBinary(), clonotype);
        }

        Map<Clonotype, Integer> samplerMap =
                sample.getClonotypes().parallelStream().collect(
                        Collectors.toMap(
                                it -> it,
                                it -> computeCorrectedCount(it, sequenceTreeMap)
                        )
                );

        return new Sample(sample, samplerMap);
    }

    private int computeCorrectedCount(Clonotype it, SequenceTreeMap<NucleotideSequence, Clonotype> sequenceTreeMap) {
        NeighborhoodIterator<NucleotideSequence, Clonotype> neighborhoodIterator =
                sequenceTreeMap.getNeighborhoodIterator(it.getCdr3ntBinary(), treeSearchParameters);

        int totalCount = it.getCount();

        for (Clonotype other : neighborhoodIterator.it()) {
            if (!other.equals(it)) {
                int mismatches = neighborhoodIterator.getMismatches();
                double logRatio = Math.log10(it.getCount() / (double) other.getCount());

                if (logRatio > mismatches * logRatioThreshold) {
                    totalCount += other.getCount();
                } else if (logRatio < -mismatches * logRatioThreshold) {

                    return 0;
                }
            }
        }

        return totalCount;

    }

    public TreeSearchParameters getTreeSearchParameters() {
        return treeSearchParameters;
    }

    public double getRatioThreshold() {
        return Math.pow(10, -logRatioThreshold);
    }
}
