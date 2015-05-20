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

import static com.antigenomics.vdjtools.Software.*
import static com.antigenomics.vdjtools.TestUtil.getResource
import static com.antigenomics.vdjtools.io.SampleStreamConnection.load

class SampleLoadTest {
    private static void loadTest(Software software, int count, int diversity) {
        def resStream = getResource("samples/${software.toString().toLowerCase()}.txt.gz")
        def sample = load(resStream, software)

        assert sample.count == count
        assert sample.diversity == diversity
    }

    @Test
    public void mitcrTest() {
        loadTest(MiTcr, 10000, 6493)
    }

    @Test
    public void migecTest() {
        loadTest(MiGec, 6147, 2420)
    }

    @Test
    public void igblastTest() {
        loadTest(IgBlast, 8408, 7274)
    }

    @Test
    public void simpleTest() {
        loadTest(VDJtools, 11878, 2257)
    }

    @Test
    public void immunoseqTest() {
        loadTest(ImmunoSeq, 10000, 2057)
    }

    @Test
    public void imgtTest() {
        loadTest(ImgtHighVQuest, 9681, 7199)
    }

    @Test
    public void mixcrTest() {
        loadTest(MiXcr, 96132, 262)
    }
}
