/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 */

package com.antigenomics.vdjtools.pool;

import com.antigenomics.vdjtools.ClonotypeWrapper;
import com.antigenomics.vdjtools.join.ClonotypeKeyGen;
import com.antigenomics.vdjtools.join.key.ClonotypeKey;
import com.antigenomics.vdjtools.sample.Clonotype;

import java.util.HashSet;
import java.util.Set;

public class StoringClonotypeAggregator extends MaxClonotypeAggregator implements ClonotypeWrapper {
    private Clonotype clonotype;
    private int convergence; 
    private final Set<ClonotypeKey> variants = new HashSet<>();
    private final ClonotypeKeyGen clonotypeKeyGen = new ClonotypeKeyGen();

    public StoringClonotypeAggregator(Clonotype clonotype, int sampleId) {
        super(clonotype, sampleId);
        this.clonotype = clonotype;
        this.convergence = 1;
        this.variants.add(clonotypeKeyGen.generateKey(clonotype));
    }

    @Override
    protected boolean tryReplace(Clonotype clonotype, int sampleId) {
        ClonotypeKey clonotypeKey = clonotypeKeyGen.generateKey(clonotype);
        
        if (!variants.contains(clonotypeKey)) {
            convergence++;
            variants.add(clonotypeKey);
        }
        
        if (super.tryReplace(clonotype, sampleId)) {
            this.clonotype = clonotype;
            return true;
        }
        return false;
    }

    public Clonotype getClonotype() {
        return clonotype;
    }

    public int getConvergence() {
        return convergence;
    }
}
