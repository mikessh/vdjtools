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
 * A class holding metadata column values for a given sample, i.e. a metadata table row 
 */
public class SampleMetadata {
    private final String sampleId
    private final ArrayList<MetadataEntry> entries
    private final MetadataTable parent

    /**
     * Creates a blank sample metadata 
     * @param sampleId unique sample id
     * @param parent parent metadata table
     */
    private SampleMetadata(String sampleId, MetadataTable parent) {
        this.sampleId = sampleId
        this.entries = new ArrayList<>()
        this.parent = parent
    }

    /**
     * Creates sample metadata and fills it with metadata entries
     * @param sampleId unique sample id
     * @param entries entries that should be assigned to a given sample
     * @param parent parent metadata table
     */
    private SampleMetadata(String sampleId, ArrayList<MetadataEntry> entries, MetadataTable parent) {
        this.sampleId = sampleId
        this.entries = entries
        this.parent = parent
    }

    /**
     * Gets the metadata entry associated with current sample and a given column
     * @param columnId metadata column id
     * @return metadata entry for a given column id
     */
    public MetadataEntry getAt(String columnId) {
        if (columnId == MetadataTable.SAMPLE_ID_COLUMN)
            return new MetadataEntry(parent, this, MetadataTable.SAMPLE_ID_COLUMN, sampleId)
        parent[sampleId, columnId]
    }

    /**
     * INTERNAL re-assigns parent metadata table 
     * @param parent
     * @return
     */
    @PackageScope
    SampleMetadata changeParent(MetadataTable parent) {
        def sampleMetadata = new SampleMetadata(sampleId, parent)
        this.entries.each { sampleMetadata.addEntry(it.changeParent(parent, sampleMetadata)) }
        sampleMetadata
    }

    /**
     * INTERNAL
     * @param sampleId
     * @param entries
     * @param parent
     * @return
     */
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

    /**
     * Gets the unique id of current sample 
     * @return sample id string
     */
    public String getSampleId() {
        sampleId
    }

    /**
     * Gets the metadata row associated with current sample
     * @return a list of metadata entries for current sample that correspond to different metadata table columns
     */
    public List<MetadataEntry> getEntries() {
        Collections.unmodifiableList(entries)
    }

    /**
     * Gets parent metadata table
     * @return metadata table this sample metadata is assigned to
     */
    public MetadataTable getParent() {
        parent
    }

    @Override
    public boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        SampleMetadata that = (SampleMetadata) o

        if (sampleId != that.sampleId) return false

        return true
    }

    @Override
    public int hashCode() {
        return sampleId.hashCode()
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        entries.size() > 0 ? entries.join("\t") : "."
    }
}
