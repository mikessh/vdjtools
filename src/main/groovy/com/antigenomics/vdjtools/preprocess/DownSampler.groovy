/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.preprocess

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.misc.MathUtil

/**
 * A class that implements down-sampling procedure, i.e.
 * selecting {@code n < N} reads from a given sample with {@code N} reads 
 */
public class DownSampler implements Sampler{
    private final Clonotype[] flattenedClonotypes
    private final Sample sample
    private final boolean unweighted

    /**
     * Create a down-sampler for the specified sample 
     * @param sample sample that would be down-sampled
     */
    public DownSampler(Sample sample) {
        this(sample, false)
    }

    /**
     * Create a down-sampler for the specified sample 
     * @param sample sample that would be down-sampled
     * @param unweighted don't weight clonotypes by frequency during sampling 
     */
    public DownSampler(Sample sample, boolean unweighted) {
        if (!unweighted && sample.count > Integer.MAX_VALUE) {
            throw new RuntimeException("Couldn't downsample samples with > ${Integer.MAX_VALUE} cells")
        }

        this.sample = sample
        this.flattenedClonotypes = new Clonotype[unweighted ? sample.diversity : sample.count]
        this.unweighted = unweighted

        int counter = 0
        sample.each {
            if (unweighted) {
                flattenedClonotypes[counter++] = it
            } else {
                for (int i = 0; i < it.count; i++) {
                    flattenedClonotypes[counter++] = it
                }
            }
        }
    }

    /**
     * Gets a specified number of reads from a given sample
     * @param count number of reads (weighted) or clonotypes (unweighted) to take
     * @return a newly create down-sampled sample, or the underlying sample if the number of reads is greated or equal to the sample size
     */
    public Sample reSample(int count) {
        if (unweighted ? count >= sample.diversity : count >= sample.count) {
            return new Sample(sample)
        } else {
            MathUtil.shuffle(flattenedClonotypes)

            def countMap = new HashMap<Clonotype, Integer>() // same as with strict overlap

            if (unweighted) {
                for (int i = 0; i < count; i++) {
                    def clonotype = flattenedClonotypes[i]
                    countMap.put(clonotype, clonotype.count)
                }
            } else {
                for (int i = 0; i < count; i++) {
                    def clonotype = flattenedClonotypes[i]
                    countMap.put(clonotype, (countMap[clonotype] ?: 0) + 1)
                }
            }

            return new Sample(sample, countMap)
        }
    }
}
