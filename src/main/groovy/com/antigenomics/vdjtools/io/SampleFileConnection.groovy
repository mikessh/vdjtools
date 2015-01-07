/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 * Last modified on 15.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata

/**
 A semi-internal class to provide lazy-loading support for SampleCollection
 */
class SampleFileConnection extends SampleStreamConnection {
    public static Sample load(String fileName, Software software) {
        load(new FileInputStreamFactory(fileName), software)
    }

    public static Sample load(String fileName, Software software, SampleMetadata sampleMetadata) {
        load(new FileInputStreamFactory(fileName), software, sampleMetadata)
    }

    public SampleFileConnection(String fileName, Software software) {
        super(new FileInputStreamFactory(fileName), software)
    }

    public SampleFileConnection(String fileName, Software software, SampleMetadata sampleMetadata) {
        this(fileName, software, sampleMetadata, false, true)
    }

    public SampleFileConnection(String fileName, Software software, SampleMetadata sampleMetadata,
                                boolean lazy, boolean store) {
        super(new FileInputStreamFactory(fileName), software, sampleMetadata, lazy, store)
    }
}
