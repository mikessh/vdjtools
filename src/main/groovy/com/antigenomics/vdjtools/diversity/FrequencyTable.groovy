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
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.join.ClonotypeKeyGen
import com.antigenomics.vdjtools.join.key.ClonotypeKey
import com.antigenomics.vdjtools.sample.Sample

public class FrequencyTable {
    private final long min, max
    private final Map<Long, Long> frequencyMap = new HashMap<>()
    private final int diversity

    FrequencyTable(Sample sample, IntersectionType intersectionType) {
        Iterable<Countable> counters

        // collapse clonotypes by a specific key
        def clonotypeKeyGen = new ClonotypeKeyGen(intersectionType)

        def hashedCounts = new HashMap<ClonotypeKey, Counter>()

        sample.each {
            def key = clonotypeKeyGen.generateKey(it)
            def counter = hashedCounts[key]
            if (!counter)
                hashedCounts.put(key, counter = new Counter())
            counter.add(it)
        }

        this.diversity = hashedCounts.size()
        counters = hashedCounts.values()

        // compute frequency table
        long min = Long.MAX_VALUE, max = -1

        counters.each {
            long count = it.count
            frequencyMap.put(count, (frequencyMap[count] ?: 0L) + 1L)
            min = Math.min(count, min)
            max = Math.max(count, max)
        }

        this.min = min
        this.max = max
    }

    public int getDiversity() {
        return diversity
    }

    public long getAt(long clonotypeSize) {
        frequencyMap[clonotypeSize] ?: 0
    }

    long getMin() {
        return min
    }

    long getMax() {
        return max
    }

    public final List<BinInfo> getBins() {
        frequencyMap.sort { it.key }.collect { new BinInfo(it.key, it.value) }
    }

    //@Override
    //public Iterator<BinInfo> iterator() {
    //    def innerIter = frequencyMap.entrySet().iterator()
    //    [hasNext: { innerIter.hasNext() },
    //     next   : { def entry = innerIter.next(); new BinInfo(entry.key, entry.value) }] as Iterator
    //}

    public static class BinInfo {
        private final long clonotypeSize, numberOfClonotypes

        BinInfo(long clonotypeSize, long numberOfClonotypes) {
            this.clonotypeSize = clonotypeSize
            this.numberOfClonotypes = numberOfClonotypes
        }

        long getClonotypeSize() {
            clonotypeSize
        }

        long getNumberOfClonotypes() {
            numberOfClonotypes
        }
    }
}
