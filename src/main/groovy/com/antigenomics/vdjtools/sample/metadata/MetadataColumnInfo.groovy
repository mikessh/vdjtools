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
 *
 * Last modified on 21.1.2015 by mikesh
 */



package com.antigenomics.vdjtools.sample.metadata

/**
 * An object containing general information for a metadata column
 */
public class MetadataColumnInfo {
    private final MetadataColumnType metadataColumnType
    private final int numericSamples, factorSamples
    private final Set<MetadataEntry> uniqueEntries = new HashSet<>()
    private final String columnId
    private final MetadataTable parent

    /**
     * Generates an info for a specified column of metadata table
     * @param table metadata table
     * @param columnId id of column to summarize
     */
    public MetadataColumnInfo(MetadataTable table, String columnId) {
        this.parent = table
        this.columnId = columnId
        int numericSamples = 0, factorSamples = 0
        table.getColumn(columnId).each {
            if (it.isNumeric())
                numericSamples++
            else
                factorSamples++
            this.uniqueEntries.add(it) // getter/member mess in groovy
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

    /**
     * Gets the type of metadata column 
     * @return numeric , semi-numeric or factor, depending on column content
     */
    public MetadataColumnType getMetadataColumnType() {
        metadataColumnType
    }

    /**
     * Gets the number of numeric samples 
     * @return number of samples that have a numeric value in a given metadata column
     */
    public int getNumericSamples() {
        numericSamples
    }

    /**
     * Gets the number of factor samples 
     * @return number of samples that have a non-numeric value in a given metadata column
     */
    public int getFactorSamples() {
        factorSamples
    }

    /**
     * Gets the set of unique entries 
     * @return and unmodifiable set contaiting unique column values 
     */
    public Set<MetadataEntry> getUniqueEntries() {
        Collections.unmodifiableSet(uniqueEntries)
    }

    /**
     * Gets the column id 
     * @return id of metadata column
     */
    public String getColumnId() {
        columnId
    }

    /**
     * Gets parent table
     * @return parent metadata table
     */
    public MetadataTable getParent() {
        parent
    }

    @Override
    public String toString() {
        "Info($columnId):\n" +
                "columnType=$metadataColumnType;\n" +
                "numericSamples=$numericSamples,factorSamples=$factorSamples;\n" +
                "uniqueEntries=${uniqueEntries.join(",")}"
    }
}
