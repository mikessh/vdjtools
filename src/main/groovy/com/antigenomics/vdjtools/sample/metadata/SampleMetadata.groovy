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
 * Last modified on 15.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.sample.metadata

import groovy.transform.PackageScope

class SampleMetadata {
    final String sampleId
    final ArrayList<MetadataEntry> entries
    final MetadataTable parent

    private SampleMetadata(String sampleId, MetadataTable parent) {
        this.sampleId = sampleId
        this.entries = new ArrayList<>()
        this.parent = parent
    }

    private SampleMetadata(String sampleId, ArrayList<MetadataEntry> entries, MetadataTable parent) {
        this.sampleId = sampleId
        this.entries = entries
        this.parent = parent
    }

    public MetadataEntry getAt(String columnId) {
        parent[sampleId, columnId]
    }

    @PackageScope
    SampleMetadata changeParent(MetadataTable parent) {
        def sampleMetadata = new SampleMetadata(sampleId, parent)
        this.entries.each { sampleMetadata.addEntry(it.changeParent(parent, sampleMetadata)) }
        sampleMetadata
    }

    @PackageScope
    static SampleMetadata predefined(String sampleId, ArrayList<MetadataEntry> entries, MetadataTable parent) {
        new SampleMetadata(sampleId, entries, parent)
    }

    @PackageScope
    void addEntry(String columnId, String value) {
        entries.add(new MetadataEntry(parent, this, columnId, value))
    }

    @PackageScope
    void addEntry(MetadataEntry metadataEntry) {
        addEntry(metadataEntry.columnId, metadataEntry.value)
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        SampleMetadata that = (SampleMetadata) o

        if (sampleId != that.sampleId) return false

        return true
    }

    @Override
    int hashCode() {
        return sampleId.hashCode()
    }

    @Override
    String toString() {
        entries.size() > 0 ? entries.join("\t") : "."
    }
}
