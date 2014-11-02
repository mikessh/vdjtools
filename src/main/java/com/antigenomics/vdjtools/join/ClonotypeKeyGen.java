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
 * Last modified on 26.10.2014 by mikesh
 */

package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.Clonotype;
import com.antigenomics.vdjtools.intersection.IntersectionType;
import com.antigenomics.vdjtools.join.key.*;
import com.antigenomics.vdjtools.sample.Sample;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.HashSet;
import java.util.Set;

public class ClonotypeKeyGen {
    private final IntersectionType intersectionType;

    public ClonotypeKeyGen(IntersectionType intersectionType) {
        this.intersectionType = intersectionType;
    }

    public Set<ClonotypeKey> generateKeySet(Sample sample) {
        Set<ClonotypeKey> keySet = new HashSet<>();
        for (Clonotype clonotype : sample) {
            keySet.add(generateKey(clonotype));
        }
        return keySet;
    }

    public ClonotypeKey generateKey(Clonotype clonotype) {
        switch (intersectionType) {
            case Nucleotide:
                return new NtKey(clonotype);

            case NucleotideV:
                return new NtVKey(clonotype);

            case NucleotideVJ:
                return new NtVJKey(clonotype);

            case AminoAcid:
                return new AaKey(clonotype);

            case AminoAcidV:
                return new AaVKey(clonotype);

            case AminoAcidVJ:
                return new AaVJKey(clonotype);

            case AminoAcidNonNucleotide:
                return new AaNotNtKey(clonotype);

            case Strict:
                return new StrictKey(clonotype);

            default:
                throw new NotImplementedException();
        }
    }
}
