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

    public boolean pass(Clonotype clonotype) {
        boolean pass = checkPass(clonotype);

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

    public int getPassedClonotypes() {
        return passedClonotypes.get();
    }

    public int getTotalClonotypes() {
        return totalClonotypes.get();
    }

    public long getPassedCount() {
        return passedCount.get();
    }

    public long getTotalCount() {
        return totalCount.get();
    }

    public double getPassedFreq() {
        return passedFreq.get();
    }

    public double getTotalFreq() {
        return totalFreq.get();
    }
}
