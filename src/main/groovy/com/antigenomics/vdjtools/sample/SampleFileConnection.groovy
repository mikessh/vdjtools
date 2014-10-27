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

package com.antigenomics.vdjtools.sample

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.metadata.SampleMetadata
import com.antigenomics.vdjtools.util.CommonUtil

/**
 A semi-internal class to provide lazy-loading support for SampleCollection
 */
class SampleFileConnection implements SampleConnection {
    private final String fileName
    private final SampleMetadata sampleMetadata
    private final Software software
    private Sample _sample = null
    private final boolean lazy, store

    SampleFileConnection(String fileName, SampleMetadata sampleMetadata, Software software,
                           boolean lazy, boolean store) {
        this.fileName = fileName
        // todo: make access to metadata as explicit as possible
        this.sampleMetadata = sampleMetadata
        this.software = software
        this.lazy = lazy
        this.store = store

        if (!lazy)
            _sample = load()
    }

    private Sample load() {
        println "[${new Date()} SampleStreamConnection] Loading sample $sampleMetadata.sampleId"
        def inputStream = CommonUtil.getFileStream(fileName)
        def sample = Sample.fromInputStream(inputStream, sampleMetadata, software)
        println "[${new Date()} SampleStreamConnection] Loaded sample $sampleMetadata.sampleId with " +
                "$sample.diversity clonotypes and $sample.count cells. " + CommonUtil.memoryFootprint()
        sample
    }

    public Sample getSample() {
        _sample ?: (store ? (_sample = load()) : load())
    }
}
