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

class AminoAcidPwm {
    private final AtomicDoubleArray pwm
    private final AtomicDouble totalFreq = new AtomicDouble();
    private final int N, M

    public AminoAcidPwm(int length) {
        this.N = length
        this.M = CommonUtil.AAS.length
        this.pwm = new AtomicDoubleArray(N * M)
    }

    private int index(int pos, byte aa) {
        if (pos < 0 || pos >= N || aa < 0 || aa >= M)
            throw new IndexOutOfBoundsException()
        aa * N + pos
    }

    public void append(Clonotype clonotype) {
        double freq = clonotype.freq
        totalFreq.addAndGet(freq)
        clonotype.cdr3aa.eachWithIndex { char aa, int pos ->
            pwm.addAndGet(index(pos, CommonUtil.aa2code(aa)), freq)
        }
    }

    public double getAt(int pos, char aa) {
        getAt(pos, CommonUtil.aa2code(aa)) / totalFreq.get()
    }

    public double getAt(int pos, byte aaCode) {
        pwm.get(index(pos, aaCode)) / totalFreq.get()
    }

    public int getN() {
        return N
    }
}
