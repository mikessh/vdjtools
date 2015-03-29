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

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.ClonotypeWrapperContainer;

import java.util.*;

public class PooledSample implements ClonotypeWrapperContainer<StoringClonotypeAggregator> {
    private final List<StoringClonotypeAggregator> clonotypes;
    private final long count;

    public PooledSample(SampleAggregator<StoringClonotypeAggregator> sampleAggregator) {
        this.clonotypes = new ArrayList<>(sampleAggregator.getDiversity());

        long count = 0;

        for (StoringClonotypeAggregator clonotypeAggregator : sampleAggregator) {
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
