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
 * Last modified on 2.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import org.apache.commons.math3.stat.regression.SimpleRegression

class FrequencyStat {
    private final List<FrequencyTable.BinInfo> binInfoList
    private final double beta, betaConf, alpha, r

    public FrequencyStat(FrequencyTable frequencyTable) {
        this.binInfoList = frequencyTable.bins
        def lm = new SimpleRegression()
        binInfoList.each {
            if (it.complementaryCdf > 0)
                lm.addData(Math.log10(it.clonotypeFreq), Math.log10(it.complementaryCdf))
        }

        this.alpha = Math.pow(10, lm.intercept)
        this.beta = lm.slope
        this.betaConf = lm.slopeConfidenceInterval
        this.r = lm.r
    }

    double getBeta() {
        beta
    }

    double getBetaConf() {
        betaConf
    }

    double getAlpha() {
        alpha
    }

    double getR() {
        r
    }
}
