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

package com.antigenomics.vdjtools.operate

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.MathUtil

/**
 * A class that implements down-sampling procedure, i.e.
 * selecting {@code n < N} reads from a given sample with {@code N} reads 
 */
public class DownSampler {
    private final Clonotype[] flattenedClonotypes
    private final Sample sample

    /**
     * Create a down-sampler for the specified sample 
     * @param sample sample that would be down-sampled
     */
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

    /**
     * Gets a specified number of reads from a given sample
     * @param count number of reads to take
     * @return a newly create down-sampled sample, or the underlying sample if the number of reads is greated or equal to the sample size
     */
    public Sample reSample(int count) {
        if (count >= sample.count) {
            return new Sample(sample)
        } else {
            MathUtil.shuffle(flattenedClonotypes)

            def countMap = new HashMap<Clonotype, Integer>() // same as with strict overlap

            for (int i = 0; i < count; i++) {
                def clonotype = flattenedClonotypes[i]
                countMap.put(clonotype, (countMap[clonotype] ?: 0) + 1)
            }

            return new Sample(sample, countMap)
        }
    }
}
