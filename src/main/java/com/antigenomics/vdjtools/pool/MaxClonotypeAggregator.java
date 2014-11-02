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

public class MaxClonotypeAggregator implements ClonotypeAggregator {
    private double maxFreq;

    public MaxClonotypeAggregator() {
        this(0.0);
    }

    public MaxClonotypeAggregator(double maxFreq) {
        this.maxFreq = maxFreq;
    }

    @Override
    public void combine(ClonotypeAggregator other) {
        this.maxFreq = Math.max(((MaxClonotypeAggregator) other).maxFreq, maxFreq);
    }

    @Override
    public void combine(Clonotype other) {
        this.maxFreq = Math.max(other.getFreq(), maxFreq);
    }

    public double getMaxFreq() {
        return maxFreq;
    }
}
