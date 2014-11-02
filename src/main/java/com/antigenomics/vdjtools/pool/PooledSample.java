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

import java.util.HashMap;
import java.util.Map;

public class PooledSample<T extends ClonotypeAggregator> {
    private final Map<ClonotypeKey, T> innerMap = new HashMap<>();
    private final ClonotypeKeyGen clonotypeKeyGen;


    public PooledSample(Sample[] samples,
                        ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory) {
        this(samples, clonotypeAggregatorFactory, IntersectionType.Strict);
    }

    public PooledSample(Sample[] samples,
                        ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory,
                        IntersectionType intersectionType) {
        this.clonotypeKeyGen = new ClonotypeKeyGen(intersectionType);
        for (Sample sample : samples) {
            for (Clonotype clonotype : sample) {
                ClonotypeKey clonotypeKey = clonotypeKeyGen.generateKey(clonotype);

                ClonotypeAggregator clonotypeAggregator = getAt(clonotypeKey);
                if (clonotypeAggregator == null) {
                    innerMap.put(clonotypeKey, clonotypeAggregatorFactory.create(clonotype));
                } else {
                    clonotypeAggregator.combine(clonotype);
                }
            }
        }
    }

    private T getAt(ClonotypeKey clonotypeKey) {
        return innerMap.get(clonotypeKey);
    }

    public T getAt(Clonotype clonotype) {
        return getAt(clonotypeKeyGen.generateKey(clonotype));
    }
}