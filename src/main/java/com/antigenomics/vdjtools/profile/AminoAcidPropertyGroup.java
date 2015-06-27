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

package com.antigenomics.vdjtools.profile;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import static com.milaboratory.core.sequence.AminoAcidSequence.ALPHABET;

public class AminoAcidPropertyGroup {
    private static final String MISSING_PROPERTY = "NA";

    private final String name;
    private final Set<String> properties = new HashSet<>();
    private final String[] propertyByAA = new String[ALPHABET.size()];

    public AminoAcidPropertyGroup(String name) {
        this.name = name;
    }

    public void put(char aminoAcid, String property) {
        properties.add(property);
        propertyByAA[ALPHABET.codeFromSymbol(aminoAcid)] = property;
    }

    public String getAt(char aminoAcid) {
        return getAt(ALPHABET.codeFromSymbol(aminoAcid));
    }

    public String getAt(byte code) {
        if (code < 0 || code >= propertyByAA.length) {
            throw new IndexOutOfBoundsException("Unknown amino acid code " + code);
        }

        String property = propertyByAA[code];
        return property.length() > 0 ? property : MISSING_PROPERTY;
    }

    public String getName() {
        return name;
    }

    public Set<String> getProperties() {
        return Collections.unmodifiableSet(properties);
    }
}
