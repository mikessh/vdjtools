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
import com.antigenomics.vdjtools.util.MathUtil

public class DownSampler {
    private final Clonotype[] flattenedClonotypes
    private final Sample sample

    public DownSampler(Sample sample) {
        if (sample.count > Integer.MAX_VALUE)
            throw new RuntimeException("Couldn't downsample samples with > ${Integer.MAX_VALUE} cells")

        this.sample = sample
        this.flattenedClonotypes = new Clonotype[sample.count]

        int counter = 0
        sample.each {
            for (int i = 0; i < it.count; i++)
                flattenedClonotypes[counter++] = it
        }
    }

    public Sample reSample(int count) {
        if (count >= sample.count) {
            return new Sample(sample)
        } else {
            MathUtil.shuffle(flattenedClonotypes)

            def countMap = new HashMap<Clonotype, Integer>() // same as with strict intersection

            for (int i = 0; i < count; i++) {
                def clonotype = flattenedClonotypes[i]
                countMap.put(clonotype, (countMap[clonotype] ?: 0) + 1)
            }

            return new Sample(sample, countMap)
        }
    }
}
