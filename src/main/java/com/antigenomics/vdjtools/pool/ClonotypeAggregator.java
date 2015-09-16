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

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.Misc;
import com.antigenomics.vdjtools.sample.Clonotype;

public abstract class ClonotypeAggregator {
    private int sampleId;
    private int incidenceCount, count;
    private double freq, freqRem = 0;

    public ClonotypeAggregator(Clonotype clonotype, int sampleId) {
        this.incidenceCount = 1;
        this.freqRem = clonotype.getFreq();
        this.count = clonotype.getCount();
        this.sampleId = sampleId;
    }

    final void combine(Clonotype other, int sampleId) {
        if (this.sampleId != sampleId) {
            this.sampleId = sampleId;
            incidenceCount++;
            freq += freqRem > 0 ? Math.log10(freqRem) : 0;  // add remainder
            freqRem = other.getFreq();                      // start accumulating freq
        } else {
            freqRem += other.getFreq();                     // continue accumulating freq
        }

        count += other.getCount();

        tryReplace(other, sampleId);
    }

    protected abstract boolean tryReplace(Clonotype other, int sampleId);

    public int getIncidenceCount() {
        return incidenceCount;
    }

    public int getCount() {
        return count;
    }

    public double getFreqLogSum() {
        return freq + (freqRem > 0 ? Math.log10(freqRem) : 0);
    }

    public double getFreqGeomMean() {
        return Math.pow(10.0, getFreqLogSum() / incidenceCount);
    }

    public double getFreqGeomMean(int numberOfSamples) {
        return Math.pow(10.0, (getFreqLogSum() + (numberOfSamples - incidenceCount) * Misc.JITTER_LOG10) / numberOfSamples);
    }
}
