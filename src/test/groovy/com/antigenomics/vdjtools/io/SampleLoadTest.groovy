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

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.Software
import org.junit.Test

import static com.antigenomics.vdjtools.TestUtil.getResource
import static com.antigenomics.vdjtools.io.SampleStreamConnection.load

class SampleLoadTest {
    private static void loadTest(String sampleName, Software software, int count, int diversity) {
        def sample = load(getResource("samples/${sampleName}.txt.gz"), software)
        assert sample.count == count
        assert sample.diversity == diversity
    }

    @Test
    public void mitcrTest() {
        loadTest("mitcr", Software.MiTcr, 10000, 6493)
    }

    @Test
    public void migecTest() {
        loadTest("migec", Software.MiGec, 6147, 2420)
    }

    @Test
    public void igblastTest() {
        loadTest("igblast", Software.IgBlast, 8408, 7274)
    }

    @Test
    public void simpleTest() {
        loadTest("simple", Software.Simple, 11878, 2257)
    }
}
