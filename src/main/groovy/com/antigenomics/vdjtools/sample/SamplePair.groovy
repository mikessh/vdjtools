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

import com.antigenomics.vdjtools.io.DummySampleConnection
import com.antigenomics.vdjtools.io.SampleConnection

class SamplePair {
    private final SampleConnection sample1conn, sample2conn
    final int i, j

    SamplePair(SampleConnection sample1conn, SampleConnection sample2conn, int i, int j) {
        this.sample1conn = sample1conn
        this.sample2conn = sample2conn
        this.i = i
        this.j = j
    }

    SamplePair(Sample sample1, Sample sample2, int i, int j) {
        this(new DummySampleConnection(sample1), new DummySampleConnection(sample2), i, j)
    }

    SamplePair(SampleConnection sample1conn, SampleConnection sample2conn) {
        this(sample1conn, sample2conn, 0, 1)
    }

    SamplePair(Sample sample1, Sample sample2) {
        this(sample1, sample2, 0, 1)
    }

    public SamplePair getReverse() {
        new SamplePair(sample2conn, sample1conn, j, i)
    }


    Sample getAt(int index) {
        switch (index) {
            case 0:
                return sample1conn.sample
            case 1:
                return sample2conn.sample
        }
        throw new IndexOutOfBoundsException()
    }
}
