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

package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.Mutation
import com.antigenomics.vdjtools.util.CommonUtil

class MotifPwm {
    final int leftSize, rightSize
    final double[][] pwmRC, pwm
    int total = 0, rgyw = 0

    MotifPwm() {
        this(2, 2)
    }

    MotifPwm(int leftSize, int rightSize) {
        this.leftSize = leftSize
        this.rightSize = rightSize
        pwmRC = new double[4][leftSize + rightSize + 1]
        pwm = new double[4][leftSize + rightSize + 1]
    }

    void add(Mutation mutation) {
        String motif = mutation.getMotif(leftSize, rightSize),
               motifRC = CommonUtil.rc(motif)
        updatePwm(pwm, motif)
        updatePwm(pwmRC, motifRC)
        total++
        if ((motif + " " + motifRC) =~ /[AG]G[TC][AT]/)
            rgyw++
    }

    private static void updatePwm(double[][] pwm, String motif) {
        motif.toCharArray().eachWithIndex { char nt, int j ->
            int code = CommonUtil.nt2code(nt)
            if (code < 0) {
                (0..3).each { int i ->
                    pwm[i][j] += 0.25
                }
            } else {
                pwm[code][j] += 1.0
            }
        }
    }

    void addAll(Collection<Mutation> mutations) {
        mutations.each { add(it) }
    }

    @Override
    String toString() {
        "NT\t" + (-leftSize..rightSize).join("\t") + "\n" +
                (0..3).collect { "${CommonUtil.code2nt(it)}\t" + pwm[it].collect().join("\t") }.join("\n") + "\n" +
        "NT_RC\t" + (-leftSize..rightSize).join("\t") + "\n" +
                (0..3).collect { "${CommonUtil.code2nt(it)}\t" + pwmRC[it].collect().join("\t") }.join("\n") + "\n" +
                "Total\t$total\nRGYW\t$rgyw"
    }
}
