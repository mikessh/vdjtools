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

class Cdr3AminoAcidProfileBuilder {
    private final int nBins
    private final boolean weighted
    private final List<String> groups

    Cdr3AminoAcidProfileBuilder(int nBins, boolean weighted, String... groups) {
        this.nBins = nBins
        this.weighted = weighted
        this.groups = groups
    }

    AminoAcidProfile create(Iterable<Clonotype> clonotypes) {
        def profile = new AminoAcidProfile(nBins,
                BasicAminoAcidProperties.INSTANCE.getGroups(groups))

        clonotypes.each {
            profile.update(it.cdr3aaBinary, weighted ? it.count : 1)
        }

        profile
    }

    int getnBins() {
        nBins
    }

    int getWeighted() {
        weighted
    }

    List<String> getGroups() {
        Collections.unmodifiableList(groups)
    }
}
