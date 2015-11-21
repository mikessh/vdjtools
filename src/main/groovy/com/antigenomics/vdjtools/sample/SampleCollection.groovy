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

package com.antigenomics.vdjtools.sample

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.io.DummySampleConnection
import com.antigenomics.vdjtools.io.SampleConnection
import com.antigenomics.vdjtools.io.SampleFileConnection
import com.antigenomics.vdjtools.sample.metadata.BlankMetadataEntryFilter
import com.antigenomics.vdjtools.sample.metadata.MetadataEntryFilter
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.sample.metadata.MetadataUtil
import com.antigenomics.vdjtools.util.ExecUtil

/**
 * Base class used in VDJtools to store and handle collections of samples.
 * Implements methods for loading clonotype tables and sample information (metadata).
 * Note that all samples in within the same sample collection should have the same {@code Software} type.
 */
class SampleCollection implements Iterable<Sample> {
    final Map<String, SampleConnection> sampleMap = new HashMap<>()
    private final Software software
    private final boolean strict, lazy, store

    private final MetadataTable metadataTable

    /**
     * Gets a metadata table that allows querying and ordering of samples in this collection.
     * @return metadata table.
     */
    public MetadataTable getMetadataTable() {
        metadataTable
    }

    protected SampleCollection(Software software, boolean strict, boolean lazy, boolean store,
                               MetadataTable metadataTable) {
        this.software = software
        this.strict = strict
        this.lazy = lazy
        this.store = store
        this.metadataTable = metadataTable
    }

    /**
     * Selects a subset of samples from this sample collection on a specified rule.
     * This method creates a deep copy of sample collection.
     * @param entryFilter sample metadata filtering rule.
     * @param sampleIds sample IDs to keep.
     * @return a copy of sample collection containing selected samples only.
     */
    public SampleCollection select(MetadataEntryFilter metadataEntryFilter = BlankMetadataEntryFilter.INSTANCE,
                                   Set<String> sampleIds = metadataTable.sampleIds) {
        def metadataTable = metadataTable.select(metadataEntryFilter, sampleIds)

        def sampleCollection = new SampleCollection(software,
                strict, lazy, store, metadataTable)

        sampleMap.each {
            def sampleConnection = it.value

            if (sampleConnection instanceof DummySampleConnection) {
                // sample is already loaded - need to reassing sample metadata
                sampleConnection = new DummySampleConnection(
                        new Sample(sampleConnection.getSample(), metadataTable.getRow(it.key))
                )
            }

            sampleCollection.sampleMap.put(it.key, sampleConnection)
        }

        sampleCollection
    }

    /**
     * Builds a sample collection from a pre-defined list of samples.
     * Samples should belong to the same metadata table, sample order will be preserved.
     * @param samples list of samples.
     */
    public static SampleCollection fromSampleList(List<Sample> samples) {
        def originalMetadata = samples[0].sampleMetadata.parent,
        // create new metadata table
            metadataTable = new MetadataTable(samples[0].sampleMetadata.parent.columnIterator.collect().toList())

        def sampleCollection = new SampleCollection(Software.VDJtools,
                true, false, true, metadataTable)

        samples.each {
            def sampleMetadata = it.sampleMetadata
            if (sampleMetadata.parent != originalMetadata)
                throw new Exception("Only samples coming from same metadata table are allowed")
            // clone metadata to a new metadata table
            sampleMetadata.changeParent(metadataTable)
            def sampleId = sampleMetadata.sampleId
            sampleCollection.sampleMap.put(sampleId, new DummySampleConnection(it))
        }
        sampleCollection
    }

    /**
     * Builds a sample collection from a pre-defined list of sample file names.
     * Samples will be assigned to generic metadata table, sample order will be preserved.
     * @param sampleFileNames list of sample file names
     * @param store if set to true, all loaded samples will be stored in memory (only has effect if lazy is set to true)
     * @param lazy if set to true, all samples will be immediately loaded, otherwise samples will be loaded by request
     * @param strict if set to false, will ignore samples with missing files, otherwise will throw an exception in such case
     * @param sort not sort sample metadata by sample id 
     */
    public SampleCollection(List<String> sampleFileNames, Software software,
                            boolean store, boolean lazy, boolean strict, boolean sort) {
        this.software = software
        this.strict = strict
        this.lazy = lazy
        this.store = store
        this.metadataTable = new MetadataTable()
        sampleFileNames.each { String fileName ->
            if (new File(fileName).exists()) {
                def sampleId = MetadataUtil.fileName2id(fileName)
                def sampleMetadata = metadataTable.createRow(sampleId)
                sampleMap.put(sampleId, new SampleFileConnection(fileName, software, sampleMetadata, lazy, store))
            } else if (strict) {
                throw new FileNotFoundException("Missing sample file $fileName")
            } else {
                println "[${new Date()} SampleCollection] WARNING: File $fileName not found, skipping"
            }
        }

        if (sort)
            metadataTable.sort()
    }

    /**
     * Builds a sample collection from a pre-defined list of sample file names.
     * Samples will be assigned to generic metadata table, sample order will be preserved.
     * @param sampleFileNames list of sample file names.
     * @param store if set to true, all loaded samples will be stored in memory (only has effect if lazy is set to true).
     * @param lazy if set to true, all samples will be immediately loaded, otherwise samples will be loaded by request.
     */
    public SampleCollection(List<String> sampleFileNames, Software software,
                            boolean store, boolean lazy) {
        this(sampleFileNames, software, store, lazy, true, false)
    }

    /**
     * Builds a sample collection from a pre-defined list of sample file names.
     * Samples will be assigned to generic metadata table, sample order will be preserved.
     * @param sampleFileNames list of sample file names.
     * @param store if set to true, all loaded samples will be stored in memory.
     */
    public SampleCollection(List<String> sampleFileNames, Software software,
                            boolean store) {
        this(sampleFileNames, software, store, true)
    }

    /**
     * Builds a sample collection from a pre-defined list of sample file names.
     * Samples will be assigned to generic metadata table, sample order will be preserved.
     * @param sampleFileNames list of sample file names.
     */
    public SampleCollection(List<String> sampleFileNames, Software software) {
        this(sampleFileNames, software, false)
    }

    /**
     * Builds a sample collection from a pre-defined list of sample file names.
     * Samples will be assigned to generic metadata table, sample order will be preserved.
     * Samples should be in VDJtools format.
     * @param sampleFileNames list of sample file names.
     */
    public SampleCollection(List<String> sampleFileNames) {
        this(sampleFileNames, Software.VDJtools)
    }

    /**
     * Loads a sample collection using custom metadata file.
     * File should contain two columns: first with file path and second with sample id.
     * Additional columns will be stored as metadata entries.
     * First line of file should contain header that includes metadata field names.
     * Samples will be ordered as they appear in file.
     * @param sampleMetadataFileName metadata file path.
     * @param software software used to get processed samples.
     * @param store if set to true, all loaded samples will be stored in memory (only has effect if lazy is set to true).
     * @param lazy if set to true, all samples will be immediately loaded, otherwise samples will be loaded by request.
     * @param strict if set to false, will ignore samples with missing files, otherwise will throw an exception in such case.
     * @param sort whether to sort sample metadata by sample id.
     */
    public SampleCollection(String sampleMetadataFileName, Software software,
                            boolean store, boolean lazy, boolean strict, boolean sort) {
        this.software = software

        this.store = store
        this.lazy = lazy
        this.strict = strict

        MetadataTable metadataTable = null

        new File(sampleMetadataFileName).withReader { reader ->
            def headerLine = reader.readLine()

            if (!headerLine.startsWith("#")) {
                throw new Exception("Metadata header should be marked with hash prefix (#)")
            }

            def metadataColumns = headerLine.split("\t")
            metadataTable = new MetadataTable(metadataColumns.size() > 2 ? metadataColumns[2..-1] : [])

            def line, splitLine
            while ((line = reader.readLine()) != null) {
                if (!line.startsWith("#")) {
                    splitLine = line.split("\t")
                    String fileName = splitLine[0], sampleId = splitLine[1]

                    fileName = ExecUtil.absoluteSamplePath(sampleMetadataFileName, fileName)

                    def entries = splitLine.length > 2 ? splitLine[2..-1] : []

                    if (new File(fileName).exists()) {
                        def sampleMetadata = metadataTable.createRow(sampleId, entries)
                        sampleMap.put(sampleId, new SampleFileConnection(fileName, software, sampleMetadata, lazy, store))
                    } else if (strict) {
                        throw new FileNotFoundException("Missing sample file $fileName")
                    } else {
                        println "[${new Date()} SampleCollection] WARNING: File $fileName not found, skipping"
                    }
                }
            }
        }

        if (sort)
            metadataTable.sort()

        this.metadataTable = metadataTable
    }

    /**
     * Loads a sample collection using custom metadata file.
     * File should contain two columns: first with file path and second with sample id
     * Additional columns will be stored as metadata entries.
     * First line of file should contain header that includes metadata field names.
     * Samples will be ordered as they appear in file
     * @param sampleMetadataFileName metadata file path
     * @param software software used to get processed samples
     * @param store if set to true, all loaded samples will be stored in memory (only has effect if lazy is set to true)
     * @param lazy if set to true, all samples will be immediately loaded, otherwise samples will be loaded by request
     */
    public SampleCollection(String sampleMetadataFileName, Software software,
                            boolean store, boolean lazy) {
        this(sampleMetadataFileName, software, store, lazy, true, false)
    }

    /**
     * Loads a sample collection using custom metadata file.
     * File should contain two columns: first with file path and second with sample id.
     * Additional columns will be stored as metadata entries.
     * First line of file should contain header that includes metadata field names.
     * Samples will be ordered as they appear in file.
     * @param sampleMetadataFileName metadata file path.
     * @param software software used to get processed samples.
     * @param store if set to true, all loaded samples will be stored in memory.
     */
    public SampleCollection(String sampleMetadataFileName, Software software,
                            boolean store) {
        this(sampleMetadataFileName, software, store, true)
    }

    /**
     * Loads a sample collection using custom metadata file.
     * File should contain two columns: first with file path and second with sample id.
     * Additional columns will be stored as metadata entries.
     * First line of file should contain header that includes metadata field names.
     * Samples will be ordered as they appear in file.
     * @param sampleMetadataFileName metadata file path.
     * @param software software used to get processed samples.
     */
    public SampleCollection(String sampleMetadataFileName, Software software) {
        this(sampleMetadataFileName, software, false)
    }

    /**
     * Loads a sample collection using custom metadata file.
     * File should contain two columns: first with file path and second with sample id
     * Additional columns will be stored as metadata entries.
     * First line of file should contain header that includes metadata field names.
     * Samples will be ordered as they appear in file.
     * Samples should be in VDJtools format.
     * @param sampleMetadataFileName metadata file path.
     */
    public SampleCollection(String sampleMetadataFileName) {
        this(sampleMetadataFileName, Software.VDJtools)
    }

    /**
     * Lists all unique sample pairs in a given collection.
     * Pairs (i, j) are chosen such as j > i, no (i, i) pairs allowed.
     * @return a list of sample pairs
     */
    public List<SamplePair> listPairs() {
        def samplePairs = new ArrayList()

        for (int i = 0; i < size(); i++)
            for (int j = i + 1; j < size(); j++)
                samplePairs.add(this[i, j])

        samplePairs
    }

    /**
     * Lists all unique sample pairs for a given sample
     * Pairs (i, j) are chosen such as j > i, no (i, i) pairs allowed
     * Sample #i will be stored into memory
     * @param i sample index.
     * @return a list of sample pairs
     */
    public List<SamplePair> listPairs(int i) {
        def samplePairs = new ArrayList<SamplePair>()

        if (i < size() - 1) {
            // Load the sample (i.e. re-wrap to dummy connection)
            def sample1conn = new DummySampleConnection(sampleMap[metadataTable.getRow(i).sampleId].sample)

            for (int j = i + 1; j < size(); j++)
                samplePairs.add(new SamplePair(sample1conn,
                        sampleMap[metadataTable.getRow(j).sampleId]))
        }

        samplePairs
    }

    /**
     * Get sample by id
     * @param sampleId sample id
     * @return sample
     */
    //Sample getAt(String sampleId) {
    //    sampleMap[sampleId]
    //}

    /**
     * Get sample by index, according to current ordering
     * @param i sample index
     * @return sample
     */
    public Sample getAt(int i) {
        if (i < 0 || i >= metadataTable.sampleCount)
            throw new IndexOutOfBoundsException()

        sampleMap[metadataTable.getRow(i).sampleId].sample
    }

    /**
     * Gets sample pair by indices, according to current ordering
     * @param i index of first sample in pair
     * @param j index of second sample in pair
     * @return sample pair
     */
    public SamplePair getAt(int i, int j) {
        new SamplePair(sampleMap[metadataTable.getRow(i).sampleId],
                sampleMap[metadataTable.getRow(j).sampleId],
                i, j)
    }

    /**
     * Gets an iterator that iterates over samples it current collection 
     * and loads them if needed
     */
    public Iterator iterator() {
        def iter = metadataTable.sampleIterator
        return [
                hasNext: {
                    iter.hasNext()
                },
                next   : {
                    def sampleId = iter.next()
                    sampleMap[sampleId].sample
                }] as Iterator
    }

    /**
     * Quickly reads all samples collecting various statistics, such as min/max read count.
     * Do not store samples in memory
     * @return
     */
    public SampleStatistics getSampleStatistics() {
        def minCount = Long.MAX_VALUE, maxCount = 0,
            minFreq = Double.MAX_VALUE, maxFreq = 0,
            minDiversity = Integer.MAX_VALUE, maxDiversity = 0
        println "[${new Date()} SampleCollection] Collecting sample statistics"
        sampleMap.each {
            def sample = it.value.haveAGlance()
            minCount = Math.min(sample.count, minCount)
            maxCount = Math.max(sample.count, maxCount)
            minFreq = Math.min(sample.freq, minFreq)
            maxFreq = Math.max(sample.freq, maxFreq)
            minDiversity = Math.min(sample.diversity, minDiversity)
            maxDiversity = Math.max(sample.diversity, maxDiversity)
        }
        new SampleStatistics(minCount, maxCount, minFreq, maxFreq, minDiversity, maxDiversity)
    }

    /**
     * Gets the number of samples
     * @return number of samples
     */
    public int size() {
        sampleMap.size()
    }

    @Override
    public String toString() {
        "samples=" + this.collect { it.sampleMetadata.sampleId }.join(",") +
                "\nmetadata=" + metadataTable.columnIterator.collect().join(",")
    }
}
