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
 * Last modified on 23.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.pwm

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.util.CommonUtil
import com.google.common.util.concurrent.AtomicDouble
import com.google.common.util.concurrent.AtomicDoubleArray

import java.util.concurrent.atomic.AtomicInteger

class CdrPwm implements Iterable<Row> {
    private final AtomicDoubleArray pwm
    private final AtomicDouble totalFreq = new AtomicDouble()
    private final AtomicInteger count = new AtomicInteger()
    private final int N, M

    public CdrPwm(int length) {
        this.N = length
        this.M = CommonUtil.AAS.length
        this.pwm = new AtomicDoubleArray(N * M)
    }

    private int index(int pos, byte aa) {
        if (pos < 0 || pos >= N || aa < 0 || aa >= M)
            throw new IndexOutOfBoundsException()
        aa * N + pos
    }

    public void update(Clonotype clonotype) {
        double freq = clonotype.freq
        totalFreq.addAndGet(freq)
        clonotype.cdr3aa.toCharArray().eachWithIndex { char aa, int pos ->
            pwm.addAndGet(index(pos, CommonUtil.aa2code(aa)), freq)
        }
        count.incrementAndGet()
    }

    public double getAt(int pos, char aa) {
        this[pos, CommonUtil.aa2code(aa)] / totalFreq.get()
    }

    public double getAt(int pos, byte aaCode) {
        pwm.get(index(pos, aaCode)) / totalFreq.get()
    }

    public MajorVariant getMajorVariant(int pos) {
        def freqs = (0..<M).collect { byte aaCode -> this[pos, aaCode] }
        def majorFreq = freqs.max(), majorCode = (byte) freqs.findIndexOf { it == majorFreq }
        new MajorVariant(pos, majorCode, majorFreq)
    }

    public CdrPattern extractConsensus(double positionalFrequencyThreshold) {
        def pattern = (0..<M).collect { int pos ->
            def major = getMajorVariant(pos)
            major.freq >= positionalFrequencyThreshold ? major.aa : "X"
        }.join("")
        new CdrPattern(pattern, getTotalFreq(), getCount())
    }

    /**
     * Gets a normalized vector of character frequencies,
     * according to WebLogo paradigm, i.e. scaled to information content
     * @param pos
     * @return vector of frequencies, ordered according to CommonUtil.AAS
     */
    public double[] getNormalizedFreqs(int pos) {
        def freqs = CommonUtil.AAS.collect { getAt(pos, it) }
        def e = 9.5 / Math.log(2) / getCount(),
            h = (double) freqs.sum { double f -> f > 0 ? -f * Math.log(f) / Math.log(2) : 0 },
            R = Math.log(20) / Math.log(2) - e - h

        freqs.collect { it * R } as double[]
    }

    public double[] getFreqs(int pos) {
        CommonUtil.AAS.collect { getAt(pos, it) }
    }

    public int getN() {
        return N
    }

    public double getTotalFreq() {
        totalFreq.get()
    }

    public int getCount() {
        count.get()
    }

    @Override
    Iterator<Row> iterator() {
        def innerIter = (0..<N).iterator()
        [hasNext: { innerIter.hasNext() },
         next   : { def pos = innerIter.next(); new Row(pos, getFreqs(pos), getNormalizedFreqs(pos)) }] as Iterator
    }

    class Row {
        final int pos
        final double[] freqs, freqsNorm

        Row(int pos, double[] freqs, double[] freqsNorm) {
            this.pos = pos
            this.freqs = freqs
            this.freqsNorm = freqsNorm
        }
    }

    class MajorVariant {
        final int pos
        final byte aaCode
        final char aa
        final double freq

        MajorVariant(int pos, byte aaCode, double freq) {
            this.pos = pos
            this.aaCode = aaCode
            this.aa = CommonUtil.code2aa(aaCode)
            this.freq = freq
        }
    }
}
