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

package com.antigenomics.vdjtools.sample

import com.antigenomics.vdjtools.io.DummySampleConnection
import com.antigenomics.vdjtools.io.SampleConnection

/**
 * A class representing a tuple of samples
 */
public class SamplePair {
    private final SampleConnection sample1conn, sample2conn
    private int i, j

    public SamplePair(SampleConnection sample1conn, SampleConnection sample2conn, int i, int j) {
        this.sample1conn = sample1conn
        this.sample2conn = sample2conn
        this.i = i
        this.j = j
    }

    /**
     * Creates a sample pair holding references to both samples and their indices in parent collection
     * @param sample1 first sample
     * @param sample2 second sample
     * @param i index of the first sample in sample collection 
     * @param j index of the second sample in sample collection
     */
    public SamplePair(Sample sample1, Sample sample2, int i, int j) {
        this(new DummySampleConnection(sample1), new DummySampleConnection(sample2), i, j)
    }

    /**
     * Creates a sample pair holding references to both samples
     * @param sample1conn an object that can be used to load the first sample
     * @param sample2conn an object that can be used to load the second sample
     */
    public SamplePair(SampleConnection sample1conn, SampleConnection sample2conn) {
        this(sample1conn, sample2conn, 0, 1)
    }

    /**
     * Creates a sample pair holding references to both samples 
     * @param sample1 first sample
     * @param sample2 second sample
     */
    public SamplePair(Sample sample1, Sample sample2) {
        this(sample1, sample2, 0, 1)
    }

    /**
     * Swaps samples
     * @return a sample pair with samples being swapped
     */
    public SamplePair getReverse() {
        new SamplePair(sample2conn, sample1conn, j, i)
    }

    /**
     * Gets the index of first sample 
     * @return the index of first sample
     */
    public getI() {
        i
    }

    /**
     * Gets the index of second sample 
     * @return the index of second sample
     */
    public getJ() {
        j
    }

    /**
     * Gets the sample that corresponds to a given index
     * @param index index of sample, {@code 0} or {@code 1}
     * @return gets the sample specified by given index
     */
    public getAt(int index) {
        switch (index) {
            case 0:
                return sample1conn.sample
            case 1:
                return sample2conn.sample
        }
        throw new IndexOutOfBoundsException()
    }

    @Override
    public String toString() {
        "SamplePair{${sample1conn},${sample2conn}}"
    }
}
