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
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.join.ClonotypeKeyGen;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.sample.Sample;

import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class SampleAggregator<T extends ClonotypeAggregator> implements Iterable<T> {
    private final Map<ClonotypeKey, T> innerMap = new HashMap<>();
    private final ClonotypeKeyGen clonotypeKeyGen;
    private final long count;

    public SampleAggregator(Iterable<Sample> samples,
                            ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory) {
        this(samples, clonotypeAggregatorFactory, OverlapType.Strict);
    }

    public SampleAggregator(Iterable<Sample> samples,
                            ClonotypeAggregatorFactory<T> clonotypeAggregatorFactory,
                            OverlapType overlapType) {
        this.clonotypeKeyGen = new ClonotypeKeyGen(overlapType);
        int sampleId = 0;
        long count = 0;

        for (Sample sample : samples) {
            System.out.println("[" + (new Date().toString()) + " " + "SamplePool] " +
                    "Pooling sample " + sample.getSampleMetadata().getSampleId());

            for (Clonotype clonotype : sample) {
                ClonotypeKey clonotypeKey = clonotypeKeyGen.generateKey(clonotype);

                ClonotypeAggregator clonotypeAggregator = getAt(clonotypeKey);
                if (clonotypeAggregator == null) {
                    innerMap.put(clonotypeKey, clonotypeAggregatorFactory.create(clonotype, sampleId));
                } else {
                    clonotypeAggregator.combine(clonotype, sampleId);
                }
            }

            count += sample.getCount( );
            sampleId++;
        }
        
        // todo: sort

        this.count = count;
    }

    private T getAt(ClonotypeKey clonotypeKey) {
        return innerMap.get(clonotypeKey);
    }

    public T getAt(Clonotype clonotype) {
        return getAt(clonotypeKeyGen.generateKey(clonotype));
    }

    public int getDiversity() {
        return innerMap.size();
    }

    public long getCount() {
        return count;
    }

    @Override
    public Iterator<T> iterator() {
        return innerMap.values().iterator();
    }
}