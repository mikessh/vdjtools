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
 * Last modified on 9.10.2014 by mikesh
 */

package com.antigenomics.vdjtools.sample.metadata

class MetadataColumnInfo {
    final MetadataColumnType metadataColumnType
    final int numericSamples, factorSamples
    final Set<MetadataEntry> uniqueEntries = new HashSet<>()
    final String columnId
    final MetadataTable parent

    MetadataColumnInfo(MetadataTable header, String columnId) {
        this.parent = header
        this.columnId = columnId
        int numericSamples = 0, factorSamples = 0
        header.getColumn(columnId).each {
            if (it.isNumeric())
                numericSamples++
            else
                factorSamples++
            uniqueEntries.add(it)
        }

        if (factorSamples > 0) {
            if (numericSamples > 0)
                metadataColumnType = MetadataColumnType.SemiNumeric
            else
                metadataColumnType = MetadataColumnType.Factor
        } else
            metadataColumnType = MetadataColumnType.Numeric

        this.numericSamples = numericSamples
        this.factorSamples = factorSamples
    }

    @Override
    public String toString() {
        "Info($columnId):\n" +
                "columnType=$metadataColumnType;\n" +
                "numericSamples=$numericSamples,factorSamples=$factorSamples;\n" +
                "uniqueEntries=${uniqueEntries.join(",")}"
    }
}
