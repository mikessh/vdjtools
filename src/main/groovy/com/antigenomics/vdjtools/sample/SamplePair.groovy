/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
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
    public Sample getI() {
        i
    }

    /**
     * Gets the index of second sample 
     * @return the index of second sample
     */
    public Sample getJ() {
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
