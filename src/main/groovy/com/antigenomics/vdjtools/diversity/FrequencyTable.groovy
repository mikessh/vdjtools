/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Countable
import com.antigenomics.vdjtools.Counter
import com.antigenomics.vdjtools.intersection.IntersectionUtil
import com.antigenomics.vdjtools.sample.Sample

class FrequencyTable {
    private final Map<Long, Long> frequencyMap = new HashMap<>()
    private final int diversity

    FrequencyTable(Sample sample, IntersectionUtil intersectionUtil) {
        Iterable<Countable> counters

        // collapse clonotypes by a specific key

        def hashedCounts = new HashMap<String, Counter>()

        sample.each {
            def key = intersectionUtil.generateKey(it)
            def counter = hashedCounts[key]
            if (!counter)
                hashedCounts.put(key, counter = new Counter())
            counter.add(it)
        }

        this.diversity = hashedCounts.size()
        counters = hashedCounts.values()

        // compute frequency table

        counters.each {
            long count = it.count
            frequencyMap.put(count, (frequencyMap[count] ?: 0L) + 1L)
        }
    }

    public int getDiversity() {
        return diversity
    }

    public long getAt(long count) {
        frequencyMap[count] ?: 0
    }
}
