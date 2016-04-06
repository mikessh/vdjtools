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

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.ClonotypeFilter;
import com.antigenomics.vdjtools.sample.Sample;

/**
 * A clonotype filter that filters out all clonotypes that do not pass a certain ratio threshold
 * when compared to the matching most abundant clonotype in another sample.
 * Used in {@link com.antigenomics.vdjtools.preprocess.Decontaminate}.
 */
public class CountFilter extends ClonotypeFilter {
    private final SampleAggregator<MaxClonotypeAggregator> sampleAggregator;
    private final double thresholdRatio;

    public CountFilter(Iterable<Sample> samples, double thresholdRatio, boolean negative) {
        super(negative);
        this.sampleAggregator = new SampleAggregator<>(samples, new MaxClonotypeAggregatorFactory());
        this.thresholdRatio = thresholdRatio;
    }

    public CountFilter(Iterable<Sample> samples, double thresholdRatio) {
        this(samples, thresholdRatio, false);
    }

    public CountFilter(Iterable<Sample> samples) {
        this(samples, 20.0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean checkPass(Clonotype clonotype) {
        MaxClonotypeAggregator aggregator = sampleAggregator.getAt(clonotype);
        return aggregator == null ||
                aggregator.getCount() < clonotype.getCount() * thresholdRatio;
    }
}
