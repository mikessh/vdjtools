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

/**
 * Base class to handle sample metadata
 */
class MetadataTable {
    /**
     * Generic metadata table. All externally created SampleMetadata objects are internally dependent on it
     */
    static final MetadataTable GENERIC_METADATA_TABLE = new MetadataTable()

    private final HashMap<String, Integer> id2index = new HashMap<>()
    private final List<String> columnIds
    private List<String> sampleOrder = new ArrayList<>()
    private final Map<String, SampleMetadata> metadataBySample = new HashMap<>()

    /**
     * Creates an empty metadata table
     */
    MetadataTable() {
        this.columnIds = new ArrayList<>()
    }

    /**
     * Creates a metadata table with pre-defined list of columns
     * @param columnIds column names
     */
    MetadataTable(List<String> columnIds) {
        this.columnIds = columnIds

        columnIds.eachWithIndex { it, ind ->
            if (id2index.containsKey(it))
                throw new IllegalArgumentException("Duplicate metadata columns not allowed")

            id2index.put(it, ind)
        }
    }

    /**
     * Reorders sample list according to provided order
     * @param sampleOrder list of sample ids in new order
     */
    void sort(List<String> sampleOrder) {
        // Perform some checks
        if (sampleOrder.size() != sampleCount)
            throw new IllegalArgumentException("Bad sample order, " +
                    "number of sample ids and sample collection size don't match")

        sampleOrder.each {
            if (!id2index.containsKey(it))
                throw new IllegalArgumentException("Bad sample order, unknown sample id $it")
        }

        this.sampleOrder = sampleOrder
    }

    /**
     * Sorts samples in metadata according to a given column.
     * Sorting is performed according to column type (Factor/Numeric/Semi-Numeric).
     * NaNs for numeric columns are always put last
     * @param columnId column identifier
     * @param descending true for descending, false for ascending sort
     */
    void sort(String columnId, boolean descending) {
        checkColumnId(columnId)

        def columnInfo = getInfo(columnId)

        switch (columnInfo.metadataColumnType) {
            case MetadataColumnType.Numeric:
            case MetadataColumnType.SemiNumeric:
                sampleOrder.sort {
                    def value = this[it, columnId].asNumeric()
                    Double.isNaN(value) ?
                            (descending ? Double.MAX_VALUE : Double.MIN_VALUE) :
                            (descending ? -value : value)
                }
                break
            case MetadataColumnType.Factor:
                sampleOrder.sort(true, new Comparator<String>() {
                    @Override
                    int compare(String o1, String o2) {
                        descending ? -o1.compareTo(o2) : o1.compareTo(o2)
                    }
                })
                break
        }
    }

    /**
     * Add metadata column
     * @param columnId name of column
     * @param values column values for all samples in collection, in corresponding order
     */
    void addColumn(String columnId, List<String> values) {
        if (values.size() != sampleCount)
            throw new Exception("Number of values should be equal to number of samples")

        addColumnId(columnId)

        sampleOrder.eachWithIndex { it, ind ->
            metadataBySample[it].entries.add(new MetadataEntry(this, columnId, values[ind]))
        }
    }

    /**
     * Add metadata column
     * @param columnId name of column
     * @param sample2value map containing column values for
     * @param strict if true all samples should be matched. Otherwise non-existing samples will be ignored,
     *               while missing samples will take "NA" value
     */
    void addColumn(String columnId, Map<String, String> sample2value, boolean strict) {
        if (strict && sampleCount != sample2value.size())
            throw new Exception("Number of values should be equal to number of samples")

        addColumnId(columnId)

        sample2value.each {
            if (metadataBySample.containsKey(it.key)) {
                metadataBySample[it.key].entries.add(new MetadataEntry(this, columnId, it.value))
            } else if (strict) {
                throw new Exception("Unknown sample $it.key")
            }
        }

        // fill missing values
        if (!strict) {
            metadataBySample.keySet().each {
                if (!sample2value.containsKey(it)) {
                    metadataBySample[it].entries.add(new MetadataEntry(this, columnId, "NA"))
                }
            }
        }
    }

    private void addColumnId(String columnId) {
        id2index.put(columnId, columnIds.size())
        columnIds.add(columnId)
    }

    private void checkColumnId(String columnId) {
        if (!containsColumn(columnId))
            throw new IllegalArgumentException("Column $columnId doesn't exist")
    }

    public int getColumnIndex(String columnId) {
        columnId = columnId.toUpperCase()
        columnIds.findIndexOf { it.toUpperCase() == columnId }
    }

    /**
     * Check if metadata table contains given column
     * @param columnId specific column
     * @return true if metadata table contains given column, otherwise false
     */
    boolean containsColumn(String columnId) {
        columnIds.contains(columnId)
    }

    /**
     * Creates a sample row, appends the metadata to table and binds sample to a given table
     * @param sampleId sample name
     * @param rowValues values of corresponding metadata fields
     * @return sample metadata instance
     */
    SampleMetadata createRow(String sampleId, List<String> rowValues) {
        // Create metadata entries, assign current metadata as parent
        def entries = new ArrayList<MetadataEntry>()

        rowValues.eachWithIndex { String value, int i ->
            entries.add(new MetadataEntry(this, columnIds[i], value))
        }

        def sampleMetadata = SampleMetadata.predefined(sampleId, entries, this)

        // Record sample to metadata if all is fine

        addSample(sampleMetadata)

        sampleMetadata
    }

    /**
     * Creates an empty sample row, appends the metadata to table and binds sample to a given table
     * @param sampleId sample name
     * @return sample metadata instance
     */
    SampleMetadata createRow(String sampleId) {
        createRow(sampleId, new ArrayList<String>())
    }

    private void addSample(SampleMetadata row) {
        // Extensive checks
        checkSample(row)

        // Record this sample

        sampleOrder.add(row.sampleId)

        metadataBySample.put(row.sampleId, row)
    }

    private void checkSample(SampleMetadata row) {
        if (metadataBySample.containsKey(row.sampleId))
            throw new IllegalArgumentException("Duplicate samples not allowed")

        if (row.entries.size() != columnIds.size())
            throw new Exception("Bad metadata row - number of entries in row " +
                    "and in parent metadata should agree")

        for (int i = 0; i < columnIds.size(); i++) {
            def rowValue = row.entries[i].columnId, metadataValue = columnIds[i]
            if (rowValue != metadataValue)
                throw new Exception("Bad metadata row - ${i}th column don't match between row " +
                        "($rowValue) and parent metadata ($metadataValue)")
        }
    }

    /**
     * Gets info on values contained in a given metadata column
     * @param columnId column name
     * @return column info
     */
    MetadataColumnInfo getInfo(String columnId) {
        new MetadataColumnInfo(this, columnId)
    }

    /**
     * Gets all metadata entries in a given column
     * @param columnId column name
     * @return metadata entries
     */
    List<MetadataEntry> getColumn(String columnId) {
        checkColumnId(columnId)

        int index = id2index[columnId]

        getColumn(index)
    }

    /**
     * Gets all metadata entries in a given column
     * @param columnIndex column index
     * @return metadata entries
     */
    List<MetadataEntry> getColumn(int columnIndex) {
        if (columnIndex < 0 || columnIndex >= columnIds.size())
            throw new IndexOutOfBoundsException()

        sampleOrder.collect {
            metadataBySample[it].entries[columnIndex]
        }
    }

    /**
     * Gets row from metadata table
     * @param sampleId sample name
     * @return metadata entries for a given sample
     */
    SampleMetadata getRow(String sampleId) {
        metadataBySample[sampleId]
    }

    /**
     * Gets row from metadata table
     * @param sampleIndex sample index
     * @return metadata entries for a given sample
     */
    SampleMetadata getRow(int sampleIndex) {
        getRow(sampleOrder[sampleIndex])
    }

    /**
     * Iterate through sample names (row names)
     * @return
     */
    Iterator<String> getSampleIterator() {
        sampleOrder.iterator()
    }

    /**
     * Iterate through metadata entry ids (column names)
     * @return
     */
    Iterator<String> getColumnIterator() {
        columnIds.iterator()
    }

    /**
     * Gets a specific metadata entry
     * @param sampleId sample name
     * @param columnId metadata id
     * @return metadata entry
     */
    MetadataEntry getAt(String sampleId, String columnId) {
        metadataBySample[sampleId].entries[id2index[columnId]]
    }

    /**
     * Gets the number of different metadata entry types (columns) in this table
     * @return number of columns in the table
     */
    int getNumberOfColumns() {
        columnIds.size()
    }

    /**
     * Gets the number of different samples (rows) in this table
     * @return number of rows in the table
     */
    int getSampleCount() {
        metadataBySample.size()
    }

    @PackageScope
    String getColumnHeader() {
        columnIds.size() > 0 ? columnIds.collect().join("\t") : "metadata_blank"
    }

    @PackageScope
    String getColumnHeader1() {
        columnIds.size() > 0 ? columnIds.collect { "1_$it" }.join("\t") : "1_metadata_blank"
    }

    @PackageScope
    String getColumnHeader2() {
        columnIds.size() > 0 ? columnIds.collect { "2_$it" }.join("\t") : "2_metadata_blank"
    }

    @Override
    String toString() {
        "\t" + columnIds.collect().join("\t") + "\n" +
                metadataBySample.collect { it.key + "\t" + it.value.entries.collect() }.join("\n")
    }
}
