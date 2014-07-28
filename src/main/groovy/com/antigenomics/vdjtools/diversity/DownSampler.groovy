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

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.sample.Sample

class DownSampler {
    final List<Clonotype> flattenedClonotypes = new ArrayList<>(1000000)
    final Sample sample

    DownSampler(Sample sample) {
        sample.clonotypes.each {
            for (long i = 0; i < it.count; i++)
                flattenedClonotypes.add(it)
        }
        this.sample = sample
    }

    Sample reSample(int count) {
        def sampledClonotypes
        if (count >= sample.count) {
            sampledClonotypes = new LinkedList<Clonotype>(sample.clonotypes)
        } else {
            Collections.shuffle(flattenedClonotypes)

            def countMap = new HashMap<Clonotype, Integer>()
            for (int i = 0; i < count; i++) {
                def clonotype = flattenedClonotypes[i]
                countMap.put(clonotype, (countMap[clonotype] ?: 0) + 1)
            }

            sampledClonotypes = new LinkedList<Clonotype>()

            countMap.each {
                sampledClonotypes.add(it.key.changeCount(it.value, count))
            }
        }

        return new Sample(sample.sampleMetadata, sampledClonotypes)
    }
}
