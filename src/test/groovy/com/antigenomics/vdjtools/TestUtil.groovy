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
 */


package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.io.InputStreamFactory
import com.antigenomics.vdjtools.sample.SampleCollection

import java.util.zip.GZIPInputStream

import static com.antigenomics.vdjtools.io.SampleStreamConnection.load

class TestUtil {
    static final SampleCollection DEFAULT_SAMPLE_COLLECTION = loadSamples(),
                                  EMPTY_SINGLE_SAMPLE = SampleCollection.fromSampleList([load(getResource("samples/empty.txt"), Software.VDJtools)])

    private static SampleCollection loadSamples() {
        def samples = Software.values().collect {
            load(getResource("samples/${it.toString().toLowerCase()}.txt.gz"), it)
        }
        samples.add(load(getResource("samples/empty.txt"), Software.VDJtools))

        SampleCollection.fromSampleList(samples)
    }


    public static InputStreamFactory getResource(String resourceName) {
        [
                create: {
                    def is = TestUtil.class.classLoader.getResourceAsStream(resourceName)
                    resourceName.endsWith(".gz") ? new GZIPInputStream(is) : is
                },
                getId : { resourceName.split("/")[-1] }
        ] as InputStreamFactory
    }
}
