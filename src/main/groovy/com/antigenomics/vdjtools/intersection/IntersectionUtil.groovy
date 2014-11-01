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
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.util.CommonUtil
import groovyx.gpars.GParsPool
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import sun.reflect.generics.reflectiveObjects.NotImplementedException

import java.util.concurrent.atomic.AtomicInteger

/**
 * A helper class for performing all sample intersection procedures in VDJtools
 */
class IntersectionUtil {
    private final int nPerms = 1000
    private final Random rnd = new Random(2106803L)
    private final IntersectionType intersectionType
    // todo: verbosity

    /**
     * Creates a helper class to perform sample intersection
     * @param intersectionType type of clonotype intersection
     */
    IntersectionUtil(IntersectionType intersectionType) {
        this.intersectionType = intersectionType
    }

    /**
     * Generates a key (e.g. TRBVx-TGTGCTAGCTGGGCATGCTTC) according to specified intersection type
     * @param clonotype clonotype to generate a key for
     * @return string key of a clonotype
     */
    String generateKey(Clonotype clonotype) {
        switch (intersectionType) {
            case IntersectionType.NucleotideV:
                return clonotype.v + "_" + clonotype.cdr3nt
            case IntersectionType.NucleotideVJ:
                return clonotype.v + "_" + clonotype.cdr3nt + "_" + clonotype.j
            case IntersectionType.AminoAcidV:
                return clonotype.v + "_" + clonotype.cdr3aa
            case IntersectionType.AminoAcidVJ:
                return clonotype.v + "_" + clonotype.cdr3aa + "_" + clonotype.j
            case IntersectionType.Nucleotide:
                return clonotype.cdr3nt
            case IntersectionType.AminoAcid:
                return clonotype.cdr3aa
            case IntersectionType.Strict:
                return clonotype.key
            // todo: special key class, non-nucleotide
            default:
                throw new NotImplementedException()
        }
    }


}
