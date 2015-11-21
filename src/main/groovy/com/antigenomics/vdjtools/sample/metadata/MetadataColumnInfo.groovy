/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */


package com.antigenomics.vdjtools.sample.metadata

/**
 * An object containing general information for a metadata column.
 */
public class MetadataColumnInfo {
    private final MetadataColumnType metadataColumnType
    private final int numericSamples, factorSamples
    private final Set<String> uniqueValues = new HashSet<>()
    private final Map<String, List<String>> sampleIdByValue = new HashMap<>()
    private final String columnId
    private final MetadataTable parent

    /**
     * Generates an info for a specified column of metadata table.
     * @param table metadata table.
     * @param columnId id of column to summarize.
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
            this.uniqueValues.add(it.value) // getter/member mess in groovy
            def sampleIds = sampleIdByValue[it.value]
            if (sampleIds == null) {
                sampleIdByValue.put(it.value, sampleIds = new ArrayList<String>())
            }
            sampleIds.add(it.parent.sampleId)
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
     * Gets the type of metadata column.
     * @return numeric , semi-numeric or factor, depending on column content.
     */
    public MetadataColumnType getMetadataColumnType() {
        metadataColumnType
    }

    /**
     * Gets the number of numeric samples.
     * @return number of samples that have a numeric value in a given metadata column.
     */
    public int getNumericSamples() {
        numericSamples
    }

    /**
     * Gets the number of factor samples.
     * @return number of samples that have a non-numeric value in a given metadata column.
     */
    public int getFactorSamples() {
        factorSamples
    }

    /**
     * Gets the set of entry values.
     * @return and unmodifiable set containing unique column values.
     */
    public Set<String> getValues() {
        Collections.unmodifiableSet(uniqueValues)
    }

    /**
     * Gets the list of sample IDs corresponding to a given entry value in this column.
     * @param value entry value to match.
     * @return list of sample identifiers.
     */
    public List<String> getSampleIds(String value) {
        Collections.unmodifiableList(sampleIdByValue[value])
    }

    /**
     * Gets the column id.
     * @return id of metadata column.
     */
    public String getColumnId() {
        columnId
    }

    /**
     * Gets parent table.
     * @return parent metadata table.
     */
    public MetadataTable getParent() {
        parent
    }

    @Override
    public String toString() {
        "Info($columnId):\n" +
                "columnType=$metadataColumnType;\n" +
                "numericSamples=$numericSamples,factorSamples=$factorSamples;\n" +
                "uniqueEntries=${_uniqueEntries.join(",")}"
    }
}
