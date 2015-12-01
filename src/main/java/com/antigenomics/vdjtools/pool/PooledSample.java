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

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.ClonotypeWrapperContainer;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.sample.SampleCollection;

import java.util.*;

public class PooledSample implements ClonotypeWrapperContainer<StoringClonotypeAggregator> {
    private final List<StoringClonotypeAggregator> clonotypes;
    private final long count;

    @SuppressWarnings("unchecked")
    public PooledSample(SampleCollection samples) {
        this(new SampleAggregator(samples,
                new StoringClonotypeAggregatorFactory(), OverlapType.Strict));
    }

    public PooledSample(SampleAggregator<StoringClonotypeAggregator> sampleAggregator) {
        this.clonotypes = new ArrayList<>(sampleAggregator.getDiversity());

        long count = 0;

        for (StoringClonotypeAggregator clonotypeAggregator : sampleAggregator) {
            clonotypeAggregator.setParent(this);
            int x = clonotypeAggregator.getCount();
            count += x;
            clonotypes.add(clonotypeAggregator);
        }

        this.count = count;

        Collections.sort(clonotypes,
                new Comparator<StoringClonotypeAggregator>() {
                    @Override
                    public int compare(StoringClonotypeAggregator o1, StoringClonotypeAggregator o2) {
                        return Integer.compare(o2.getCount(), o1.getCount()); // inverse - sort descending
                    }
                });
    }


    @Override
    public double getFreq() {
        return 1.0;
    }

    @Override
    public long getCount() {
        return count;
    }

    @Override
    public int getDiversity() {
        return clonotypes.size();
    }

    @Override
    public StoringClonotypeAggregator getAt(int index) {
        return clonotypes.get(index);
    }

    @Override
    public boolean isSorted() {
        return true;
    }

    @Override
    public Iterator<StoringClonotypeAggregator> iterator() {
        return clonotypes.iterator();
    }
}
