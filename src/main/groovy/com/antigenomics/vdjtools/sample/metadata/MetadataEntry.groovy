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

package com.antigenomics.vdjtools.sample.metadata

import groovy.transform.PackageScope

class MetadataEntry {
    public final MetadataTable grandParent
    public final SampleMetadata parent
    public final String columnId
    public String value

    MetadataEntry(MetadataTable grandParent, SampleMetadata parent, String columnId, String value) {
        this.grandParent = grandParent
        this.parent = parent
        this.columnId = columnId
        this.value = value
    }

    @PackageScope
    public MetadataEntry changeParent(MetadataTable grandParent, SampleMetadata parent) {
        new MetadataEntry(grandParent, parent, columnId, value)
    }

    boolean isNumeric() {
        value.isDouble()
    }

    double asNumeric() {
        numeric ? value.toDouble() : Double.NaN
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        MetadataEntry that = (MetadataEntry) o

        if (columnId != that.columnId) return false
        if (value != that.value) return false

        return true
    }

    @Override
    int hashCode() {
        int result
        result = columnId.hashCode()
        result = 31 * result + value.hashCode()
        return result
    }

    @Override
    String toString() {
        value
    }
}
