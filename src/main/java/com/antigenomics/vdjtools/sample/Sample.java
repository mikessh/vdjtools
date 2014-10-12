/**
 * Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 */

package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.Clonotype;
import com.antigenomics.vdjtools.Software;
import com.antigenomics.vdjtools.parser.ClonotypeStreamParser;
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata;

import java.io.InputStream;
import java.util.*;

public class Sample implements Iterable<Clonotype> {
    private final List<Clonotype> clonotypes = new ArrayList<>();
    private final SampleMetadata sampleMetadata;
    private long count = 0;

    private Sample(SampleMetadata sampleMetadata) {
        this.sampleMetadata = sampleMetadata;
    }

    public Sample(Sample toClone, HashMap<Clonotype, Integer> samplerMap) {
        this.sampleMetadata = toClone.sampleMetadata;

        for (Clonotype clonotype : toClone.clonotypes) {
            Integer newCount = samplerMap.get(clonotype);

            if (newCount != null)
                this.addClonotype(new Clonotype(clonotype, this, newCount));
        }

        Collections.sort(clonotypes);
    }

    public Sample(Sample toClone, ClonotypeFilter filter, int top) {
        this.sampleMetadata = toClone.sampleMetadata;

        for (Clonotype clonotype : toClone.clonotypes) {
            if (top > -1 && this.getDiversity() == top)
                break;

            if (filter.pass(clonotype))
                this.addClonotype(new Clonotype(clonotype, this));
        }
    }

    public Sample(Sample toClone, ClonotypeFilter filter) {
        this(toClone, filter, -1);
    }

    public Sample(Sample toClone) {
        this(toClone, BlankClonotypeFilter.INSTANCE, -1);
    }

    public static Sample fromInputStream(InputStream inputStream,
                                         SampleMetadata sampleMetadata,
                                         Software software,
                                         int top) {
        Sample sample = new Sample(sampleMetadata);

        ClonotypeStreamParser clonotypeStreamParser = ClonotypeStreamParser.create(inputStream, software, sample);

        for (Clonotype clonotype : clonotypeStreamParser) {
            if (top > -1 && sample.getDiversity() == top)
                break;

            sample.addClonotype(clonotype);
        }

        Collections.sort(sample.clonotypes);

        return sample;
    }

    public static Sample fromInputStream(InputStream inputStream,
                                         SampleMetadata sampleMetadata,
                                         Software software) {
        return fromInputStream(inputStream, sampleMetadata, software, -1);
    }

    private void addClonotype(Clonotype clonotype) {
        count += clonotype.getCount();
        clonotypes.add(clonotype);
    }

    public SampleMetadata getSampleMetadata() {
        return sampleMetadata;
    }

    //todo: rename to total

    public long getCount() {
        return count;
    }

    public long getDiversity() {
        return clonotypes.size();
    }

    @Override
    public Iterator<Clonotype> iterator() {
        return clonotypes.iterator();
    }
}
