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

package com.antigenomics.vdjtools.timecourse

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.MutationSet

/**
 * A class representing a trace of a certain clonotype in multiple samples
 */
class DynamicClonotype {
    static final double JITTER = 1e-7
    private final Clonotype[] instances
    private double[] frequencies = null
    private Clonotype representative = null

    /**
     * Creates a dynamic clonotype from an array of clonotype instances,
     * presumably the same clonotype as observed in different samples
     * @param instances an array of clonotype instances
     */
    DynamicClonotype(Clonotype[] instances) {
        this.instances = instances
    }

    /**
     * Gets a representative clonotype, based on max frequency
     * @return representative clonotype
     */
    Clonotype getRepresentative() {
        // lazy get
        representative ?: (representative = instances.max { frequency(it) })
    }

    /**
     * Gets the vector of clonotype frequencies.
     * Missing and zero-frequency clonotypes are replaced by JITTER values.
     * @return vector of clonotype frequencies
     */
    double[] getFrequencies() {
        // lazy get
        frequencies ?: (frequencies = instances.collect { frequency(it) } as double[])
    }

    /**
     * Gets the vector of mutation sets.
     * Missing and zero-frequency clonotypes are replaced by an empty sets.
     * @return vector of mutation sets
     */
    MutationSet[] getMutations() {
        instances.collect {
            present(it) ? new MutationSet(it.mutations) : new MutationSet()
        }
    }

    /**
     * Gets harmonic mean of clonotype frequency.
     * Missing and zero-frequency clonotypes are replaced by JITTER values.
     * @return mean frequency
     */
    double getMeanFrequency() {
        double meanFreq = 1
        instances.each { meanFreq *= frequency(it) }
        Math.pow(meanFreq, 1.0 / instances.size())
    }

    private static boolean present(Clonotype clonotype) {
        clonotype && clonotype.freq > 0
    }

    private static double frequency(Clonotype clonotype) {
        present(clonotype) ? clonotype.freq : JITTER
    }

    static final List<String> PRINT_FIELDS = ["cdr3aa", "cdr3nt",
                                              "v", "d", "j",
                                              "inFrame", "isComplete", "noStop"]

    @Override
    String toString() {
        PRINT_FIELDS.collect { getRepresentative()."$it" }.join("\t")
    }
}
