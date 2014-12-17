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
 * Last modified on 17.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.sample


class SampleStatistics {
    public final long minCount, maxCount
    public final double minFreq, maxFreq
    public final int minDiversity, maxDiversity

    SampleStatistics(long minCount, long maxCount,
                     double minFreq, double maxFreq,
                     int minDiversity, int maxDiversity) {
        this.minCount = minCount
        this.maxCount = maxCount
        this.minFreq = minFreq
        this.maxFreq = maxFreq
        this.minDiversity = minDiversity
        this.maxDiversity = maxDiversity
    }

    public long getMinCount() {
        minCount
    }

    public long getMaxCount() {
        maxCount
    }

    public double getMinFreq() {
        minFreq
    }

    public double getMaxFreq() {
        maxFreq
    }

    public int getMinDiversity() {
        minDiversity
    }

    public int getMaxDiversity() {
        maxDiversity
    }
}
