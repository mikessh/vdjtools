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
import com.antigenomics.vdjtools.Util

class MotifPwm {
    final int leftSize, rightSize
    final double[][] pwmC, pwmF

    MotifPwm() {
        this(3, 3)
    }

    MotifPwm(int leftSize, int rightSize) {
        this.leftSize = leftSize
        this.rightSize = rightSize
        pwmC = new double[4][leftSize + rightSize + 1]
        pwmF = new double[4][leftSize + rightSize + 1]
    }

    void add(Mutation mutation) {
        String motif = mutation.getMotif(leftSize, rightSize)
        motif.toCharArray().eachWithIndex { char nt, int j ->
            int code = Util.nt2code(nt)
            if (code < 0) {
                (0..3).each { int i ->
                    pwmC[i][j] += 0.25
                    pwmF[i][j] += 0.25 * mutation.freq
                }
            } else {
                pwmC[code][j] += 1.0
                pwmF[code][j] += mutation.freq
            }
        }
    }

    void addAll(Collection<Mutation> mutations) {
        mutations.each { add(it) }
    }

    @Override
    String toString() {
        "NT\t" + (-leftSize..rightSize).join("\t") + "\n" +
                (0..3).collect { "${Util.code2nt(it)}\t" + pwmC[it].collect().join("\t") }.join("\n")
    }
}
