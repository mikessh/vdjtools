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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.ClonotypeWrapper
import sun.reflect.generics.reflectiveObjects.NotImplementedException

class IntersectionUtil {
    private final IntersectionType intersectionType

    IntersectionUtil(IntersectionType intersectionType) {
        this.intersectionType = intersectionType
    }

    String generateKey(Clonotype clonotype) {
        switch (intersectionType) {
            case IntersectionType.NucleotideV:
                return clonotype.cdr3nt + "\t" + clonotype.v
            case IntersectionType.Nucleotide:
                return clonotype.cdr3nt
            case IntersectionType.AminoAcid:
                return clonotype.cdr3aa
            default:
                throw new NotImplementedException()
        }
    }

    ClonotypeWrapper wrap(Clonotype clonotype) {
        new Wrapper(clonotype)
    }

    private class Wrapper implements ClonotypeWrapper{
        final Clonotype clonotype

        Wrapper(Clonotype clonotype) {
            this.clonotype = clonotype
        }

        boolean equals(o) {
            if (this.is(o)) return true

            Wrapper clonotypeWrapper = (Wrapper) o

            switch (intersectionType) {
                case IntersectionType.Nucleotide:
                    if (clonotype.cdr3nt != clonotypeWrapper.clonotype.cdr3nt)
                        return false
                    break
                case IntersectionType.NucleotideV:
                    if (clonotype.cdr3nt != clonotypeWrapper.clonotype.cdr3nt ||
                            clonotype.v != clonotypeWrapper.clonotype.v)
                        return false
                    break
                case IntersectionType.AminoAcid:
                    if (clonotype.cdr3aa != clonotypeWrapper.clonotype.cdr3aa)
                        return false
                    break
                case IntersectionType.AminoAcidNonNucleotide:
                    if (clonotype.cdr3aa != clonotypeWrapper.clonotype.cdr3aa ||
                            clonotype.cdr3nt == clonotypeWrapper.clonotype.cdr3nt)
                        return false
                    break
            }

            return true
        }

        int hashCode() {
            if (intersectionType == IntersectionType.NucleotideV)
                return clonotype.cdr3nt.hashCode() + 31 * clonotype.v.hashCode()
            else
                return intersectionType == IntersectionType.Nucleotide ?
                        clonotype.cdr3nt.hashCode() : clonotype.cdr3aa.hashCode()
        }
    }
}
