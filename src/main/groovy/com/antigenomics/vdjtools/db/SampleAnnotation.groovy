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

package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.CommonUtil

class SampleAnnotation {
    private final Sample sample
    private final sampleCdrFreqs = new HashMap<String, Double>()
    private final boolean oneMM = true

    SampleAnnotation(Sample sample) {
        this.sample = sample
        sample.each {
            sampleCdrFreqs.put(it.cdr3aa, (sampleCdrFreqs[it.cdr3aa] ?: 0) + it.freq)
        }
    }

    HashMap<String, Double> getEntryFrequencies(CdrDatabase database) {
        def dbCdrFreqs = new HashMap<String, Double>()

        database.each { String seq ->
            dbCdrFreqs.put(seq, (Double) 0)

            def freq = sampleCdrFreqs[seq]
            if (freq)
                dbCdrFreqs.put(seq, dbCdrFreqs[seq] + freq)

            if (oneMM) {
                char[] aaArr = seq.toCharArray()
                for (int i = 0; i < aaArr.length; i++) {
                    char oldAA = aaArr[i]
                    CommonUtil.AAS.each { char newAA ->
                        if (newAA != oldAA) {
                            aaArr[i] = newAA
                            def newSeq = new String(aaArr)
                            freq = sampleCdrFreqs[newSeq]
                            if (freq) {
                                // todo: slow impl, for now when db is small doesn't matter
                                dbCdrFreqs.put(seq, dbCdrFreqs[seq] + freq)
                            }
                        }
                    }
                    aaArr[i] = oldAA
                }
            }
        }

        dbCdrFreqs
    }
}
