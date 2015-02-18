/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 *
 * Last modified on 20.1.2015 by mikesh
 */



package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.sample.Sample

/**
 * A wrapper for sample stored in memory
 */
public class DummySampleConnection implements SampleConnection {
    private final Sample sample

    /**
     * Wraps a given sample into a dummy sample connection
     * @param sample sample object which is assumed to be filled with clonotypes
     */
    public DummySampleConnection(Sample sample) {
        this.sample = sample
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Sample getSample() {
        sample
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Sample haveAGlance() {
        sample
    }

    @Override
    public String toString() {
        "DummyFileConnection{.>${sample.sampleMetadata.sampleId}}"
    }
}
