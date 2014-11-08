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
 * Last modified on 2.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.Clonotype;
import com.antigenomics.vdjtools.intersection.IntersectionType;
import com.antigenomics.vdjtools.join.ClonotypeKeyGen;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.sample.Sample;

import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class SampleAggregator<T extends ClonotypeAggregator> implements Iterable<T> {
    private final Map<ClonotypeKey, T> innerMap = new HashMap<>();
    private final ClonotypeKeyGen clonotypeKeyGen;
    private final long count;

    public SampleAggregator(Iterable<Sample> samples,
                            ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory) {
        this(samples, clonotypeAggregatorFactory, IntersectionType.Strict);
    }

    public SampleAggregator(Iterable<Sample> samples,
                            ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory,
                            IntersectionType intersectionType) {
        this.clonotypeKeyGen = new ClonotypeKeyGen(intersectionType);
        int sampleId = 0;
        long count = 0;
        for (Sample sample : samples) {
            System.out.println("[" + (new Date().toString()) + " " + "SamplePool] " +
                    "Pooling sample " + sample.getSampleMetadata().getSampleId());

            for (Clonotype clonotype : sample) {
                ClonotypeKey clonotypeKey = clonotypeKeyGen.generateKey(clonotype);

                ClonotypeAggregator clonotypeAggregator = getAt(clonotypeKey);
                if (clonotypeAggregator == null) {
                    innerMap.put(clonotypeKey, clonotypeAggregatorFactory.create(clonotype, sampleId));
                } else {
                    clonotypeAggregator.combine(clonotype, sampleId);
                }
            }

            count += sample.getCount( );
            sampleId++;
        }

        this.count = count;
    }

    private T getAt(ClonotypeKey clonotypeKey) {
        return innerMap.get(clonotypeKey);
    }

    public T getAt(Clonotype clonotype) {
        return getAt(clonotypeKeyGen.generateKey(clonotype));
    }

    public int getDiversity() {
        return innerMap.size();
    }

    public long getCount() {
        return count;
    }

    @Override
    public Iterator<T> iterator() {
        return innerMap.values().iterator();
    }
}