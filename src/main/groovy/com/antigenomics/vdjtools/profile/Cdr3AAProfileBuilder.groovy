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

package com.antigenomics.vdjtools.profile

import com.antigenomics.vdjtools.sample.Clonotype

class Cdr3AAProfileBuilder {
    private final Map<Cdr3Region, Integer> binning
    private final boolean weighted
    private final List<String> groups

    Cdr3AAProfileBuilder(Map<Cdr3Region, Integer> binning, boolean weighted, String... groups) {
        this.binning = binning
        this.weighted = weighted
        this.groups = groups
    }

    Map<Cdr3Region, AminoAcidProfile> create(Iterable<Clonotype> clonotypes) {
        def profiles = new HashMap<Cdr3Region, AminoAcidProfile>()

        binning.each {
            profiles.put(it.key, new AminoAcidProfile(it.value,
                    BasicAminoAcidProperties.INSTANCE.getGroups(groups)))
        }

        clonotypes.each { clonotype ->
            if (clonotype.isCoding()) {
                profiles.each {
                    def aaSeq = it.key.extractAminoAcid(clonotype)
                    if (aaSeq.size() > 0) {
                        it.value.update(aaSeq,
                                weighted ? clonotype.count : 1)
                    }
                }
            }
        }

        profiles
    }
}
