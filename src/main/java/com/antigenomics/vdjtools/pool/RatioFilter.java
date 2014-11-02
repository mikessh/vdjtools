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
import com.antigenomics.vdjtools.sample.ClonotypeFilter;
import com.antigenomics.vdjtools.sample.Sample;

public class RatioFilter extends ClonotypeFilter {
    private final PooledSample<MaxClonotypeAggregator> pooledSample;
    private final double thresholdRatio;

    public RatioFilter(Sample[] samples, double thresholdRatio) {
        this.pooledSample = new PooledSample<>(samples, new MaxClonotypeAggregatorFactory());
        this.thresholdRatio = thresholdRatio;
    }

    public RatioFilter(Sample[] samples) {
        this(samples, 20.0);
    }

    @Override
    protected boolean checkPass(Clonotype clonotype) {
        return pooledSample.getAt(clonotype).getMaxFreq() < clonotype.getFreq() * thresholdRatio;
    }
}
