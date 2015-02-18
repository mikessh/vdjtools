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

import groovy.transform.PackageScope

/**
 * Metadata entry, an object holding metadata table value that is associated with a given sample and column 
 */
public class MetadataEntry {
    public final MetadataTable grandParent
    public final SampleMetadata parent
    public final String columnId
    public String value

    /**
     * Creates new metadata entry
     * @param grandParent metadata table
     * @param parent sample metadata, a set of column values for a given sample (e.g. metadata table row)
     * @param columnId id of metadata table column
     * @param value metadata table value that is associated with a given sample and column
     */
    public MetadataEntry(MetadataTable grandParent, SampleMetadata parent, String columnId, String value) {
        this.grandParent = grandParent
        this.parent = parent
        this.columnId = columnId
        this.value = value
    }

    /**
     * INTERNAL re-assignes a given entry to new metadata table and sample 
     * @param grandParent
     * @param parent
     * @return
     */
    @PackageScope
    MetadataEntry changeParent(MetadataTable grandParent, SampleMetadata parent) {
        new MetadataEntry(grandParent, parent, columnId, value)
    }

    /**
     * Checks if a given metadata entry contains a numeric value
     * @return {@code true} if value could be converted to number, {@code false} otherwise
     */
    public boolean isNumeric() {
        value.isDouble()
    }

    /**
     * Converts the metadata entry value to double-precision number
     * @return a double precision value if sample is numeric, {@code Double.NaN} otherwise
     */
    public double asNumeric() {
        numeric ? value.toDouble() : Double.NaN
    }

    @Override
    public boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        MetadataEntry that = (MetadataEntry) o

        if (columnId != that.columnId) return false
        if (value != that.value) return false

        return true
    }

    @Override
    public int hashCode() {
        int result
        result = columnId.hashCode()
        result = 31 * result + value.hashCode()
        return result
    }

    @Override
    public String toString() {
        value
    }
}
