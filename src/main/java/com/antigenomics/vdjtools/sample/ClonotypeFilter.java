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
 * Last modified on 17.10.2014 by mikesh
 */

package com.antigenomics.vdjtools.sample;

import com.antigenomics.vdjtools.Clonotype;
import com.google.common.util.concurrent.AtomicDouble;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

public abstract class ClonotypeFilter {
    private final AtomicInteger passedClonotypes = new AtomicInteger(),
            totalClonotypes = new AtomicInteger();
    private final AtomicLong passedCount = new AtomicLong(),
            totalCount = new AtomicLong();
    private final AtomicDouble passedFreq = new AtomicDouble(),
            totalFreq = new AtomicDouble();
    private final boolean negative;

    protected ClonotypeFilter(boolean negative) {
        this.negative = negative;
    }

    public boolean pass(Clonotype clonotype) {
        boolean pass = negative ^ checkPass(clonotype);

        if (pass) {
            passedClonotypes.incrementAndGet();
            passedCount.addAndGet(clonotype.getCount());
            passedFreq.addAndGet(clonotype.getFreq());
        }

        totalClonotypes.incrementAndGet();
        totalCount.addAndGet(clonotype.getCount());
        totalFreq.addAndGet(clonotype.getFreq());

        return pass;
    }

    protected abstract boolean checkPass(Clonotype clonotype);

    public ClonotypeFilterStats getStats() {
        return new ClonotypeFilterStats(passedClonotypes.get(), totalClonotypes.get(),
                passedCount.get(), totalCount.get(),
                passedFreq.get(), totalFreq.get());
    }

    public ClonotypeFilterStats getStatsAndFlush() {
        ClonotypeFilterStats stats = getStats();
        passedClonotypes.set(0);
        totalClonotypes.set(0);
        passedCount.set(0L);
        totalCount.set(0L);
        passedFreq.set(0.0);
        totalFreq.set(0.0);
        return stats;
    }

    public static class ClonotypeFilterStats {
        private final int passedClonotypes, totalClonotypes;
        private final long passedCount, totalCount;
        private final double passedFreq, totalFreq;

        public ClonotypeFilterStats(int passedClonotypes, int totalClonotypes,
                                    long passedCount, long totalCount,
                                    double passedFreq, double totalFreq) {
            this.passedClonotypes = passedClonotypes;
            this.totalClonotypes = totalClonotypes;
            this.passedCount = passedCount;
            this.totalCount = totalCount;
            this.passedFreq = passedFreq;
            this.totalFreq = totalFreq;
        }

        public int getPassedClonotypes() {
            return passedClonotypes;
        }

        public int getTotalClonotypes() {
            return totalClonotypes;
        }

        public long getPassedCount() {
            return passedCount;
        }

        public long getTotalCount() {
            return totalCount;
        }

        public double getPassedFreq() {
            return passedFreq;
        }

        public double getTotalFreq() {
            return totalFreq;
        }
    }
}
