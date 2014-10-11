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

package com.antigenomics.vdjtools;

import com.antigenomics.vdjtools.parser.ClonotypeParser;
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class SampleJ implements Iterable<ClonotypeJ> {
    private final List<ClonotypeJ> clonotypes = new ArrayList<>();
    private final SampleMetadata sampleMetadata;
    private long totalCount = 0;

    private SampleJ(SampleMetadata sampleMetadata) {
        this.sampleMetadata = sampleMetadata;
    }

    public SampleJ(SampleJ parent, ClonotypeFilter filter, int top) {
        this.sampleMetadata = parent.sampleMetadata;

        for (ClonotypeJ clonotype : clonotypes) {
            if (top > -1 && this.getDiversity() == top)
                break;

            if (filter.pass(clonotype))
                this.addClonotype(new ClonotypeJ(clonotype, this));
        }
    }

    public SampleJ(SampleJ parent, ClonotypeFilter filter) {
        this(parent, filter, -1);
    }

    public static SampleJ fromStream(InputStream inputStream,
                                     SampleMetadata sampleMetadata,
                                     Software software,
                                     int top) {
        SampleJ sample = new SampleJ(sampleMetadata);

        ClonotypeParser clonotypeParser = ClonotypeParser.create(inputStream, software, sample);

        for (ClonotypeJ clonotype : clonotypeParser) {
            if (top > -1 && sample.getDiversity() == top)
                break;

            sample.addClonotype(clonotype);
        }

        Collections.sort(sample.clonotypes);

        return sample;
    }

    public static SampleJ fromStream(InputStream inputStream,
                                     SampleMetadata sampleMetadata,
                                     Software software) {
        return fromStream(inputStream, sampleMetadata, software, -1);
    }

    private void addClonotype(ClonotypeJ clonotype) {
        totalCount += clonotype.getCount();
        clonotypes.add(clonotype);
    }

    public SampleMetadata getSampleMetadata() {
        return sampleMetadata;
    }

    public long getTotalCount() {
        return totalCount;
    }

    public long getDiversity() {
        return clonotypes.size();
    }

    @Override
    public Iterator<ClonotypeJ> iterator() {
        return clonotypes.iterator();
    }
}
