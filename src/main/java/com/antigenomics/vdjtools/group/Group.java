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

package com.antigenomics.vdjtools.group;

import com.antigenomics.vdjtools.sample.Clonotype;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

public class Group<SignatureType extends GroupSignature> {
    private final List<Clonotype> clonotypes = new LinkedList<>();
    private final SignatureType signature;

    public Group(SignatureType signature) {
        this.signature = signature;
    }

    public void add(Clonotype clonotype) {
        clonotypes.add(clonotype);
    }

    public SignatureType getSignature() {
        return signature;
    }

    public List<Clonotype> getClonotypes() {
        return Collections.unmodifiableList(clonotypes);
    }

    public double getAbundancePValue(Clonotype clonotype) {
        double pValue = 0;
        for (Clonotype other : clonotypes) {
            if (clonotype.getFreq() < other.getFreq()) {
                pValue++;
            } else if (clonotype.getFreq() == other.getFreq()) {
                pValue += 0.5;
            }
        }
        return pValue / clonotypes.size();
    }

    public int size() {
        return clonotypes.size();
    }
}
