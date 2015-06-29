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

import com.antigenomics.vdjtools.util.CommonUtil

class BasicAminoAcidProperties {
    static final BasicAminoAcidProperties INSTANCE = new BasicAminoAcidProperties()

    private final AminoAcidProperty[] aminoAcidProperties

    private BasicAminoAcidProperties() {
        aminoAcidProperties = AminoAcidProperty.fromInput(CommonUtil.resourceStream("profile/aa_property_table.txt"))
    }

    List<String> getPropertyNames() {
        aminoAcidProperties.collect { it.name }
    }

    AminoAcidProperty[] getProperties(List<String> propertyNames) {
        if (propertyNames.isEmpty())
            return getProperties()

        aminoAcidProperties.findAll { propertyNames.contains(it.name) } as AminoAcidProperty[]
    }

    AminoAcidProperty[] getProperties() {
        Arrays.copyOf(aminoAcidProperties, aminoAcidProperties.length)
    }
}
