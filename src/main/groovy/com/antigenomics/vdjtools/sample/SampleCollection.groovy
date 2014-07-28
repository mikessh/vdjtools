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

package com.antigenomics.vdjtools.sample

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.util.CommonUtil

/**
 * Base class to store and handle collections of samples in VDJtools
 * Provides loading samples and their metadata
 * The only limitation is that samples should be processed by same software
 */
class SampleCollection implements Iterable<Sample> {
    private final Map<String, Sample> sampleMap = new HashMap<>()
    private final Map<String, List<String>> filesBySample = new HashMap<>()
    private final HashSet<String> loadedSamples = new HashSet<>()
    private final Software software
    private final boolean strict, lazy

    private final MetadataTable metadataTable

    /**
     * Gets a metadata table that allows querying and ordering of samples in this collection
     * @return metadata table
     */
    MetadataTable getMetadataTable() {
        metadataTable
    }

    /**
     * Builds a sample collection from a pre-defined list of samples.
     * Samples should belong to the same metadata table, sample order will be preserved.
     * @param samples list of samples
     */
    SampleCollection(List<Sample> samples) {
        this.software = null
        this.strict = true
        this.lazy = false
        this.metadataTable = samples[0].sampleMetadata.parent
        samples.each {
            if (it.sampleMetadata.parent != metadataTable)
                throw new Exception("Only samples coming from same metadata table are allowed")
            def sampleId = it.sampleMetadata.sampleId
            sampleMap.put(sampleId, it)
        }
    }

    /**
     * Loads a sample collection using custom metadata file.
     * File should contain two columns: first with file path and second with sample id
     * Additional columns will be stored as metadata entries.
     * First line of file should contain header that includes metadata field names.
     * Samples will be ordered as they appear in file
     * @param sampleMetadataFileName metadata file path
     * @param software software used to get processed samples
     * @param strict corresponding files should exist for all samples in metadata
     * @param lazy if true will load samples on-demand. Otherwise pre-loads all samples
     */
    SampleCollection(String sampleMetadataFileName, Software software, boolean strict, boolean lazy) {
        this.software = software
        this.strict = strict
        this.lazy = lazy

        def nSamples = 0, nClonotypes = 0

        def metadataTable

        new File(sampleMetadataFileName).withReader { reader ->
            metadataTable = new MetadataTable(reader.readLine().split("\t")[2..-1])

            def line, splitLine
            while ((line = reader.readLine()) != null) {
                splitLine = line.split("\t")
                String fileName = splitLine[0], sampleId = splitLine[1]

                def entries = splitLine.length > 2 ? splitLine[2..-1] : []

                // This one is catched later
                //if (entries.size() != metadata.size())
                //    throw new Exception("Different number of entries in metadata header and sample $sampleId")

                def inputFile = new File(fileName)

                if (inputFile.exists()) {
                    def clonotypes = new ArrayList<Clonotype>()

                    if (lazy) {
                        def fileList = filesBySample[sampleId]
                        if (!fileList)
                            filesBySample.put(sampleId, fileList = new LinkedList<String>())
                        fileList.add(fileName)
                    } else {
                        clonotypes = SampleUtil.loadClonotypes(fileName, software)
                    }

                    def sample = sampleMap[sampleId]
                    if (!sample) {
                        def sampleMetadata = metadataTable.createSample(sampleId, entries)

                        sampleMap.put(sampleId, new Sample(sampleMetadata, clonotypes))
                    } else {
                        sample.clonotypes.addAll(clonotypes)
                    }

                    nSamples++

                    if (!lazy) {
                        nClonotypes += clonotypes.size()

                        println "[${new Date()} SampleCollection] Loaded $nSamples samples and " +
                                "$nClonotypes clonotypes so far. " + CommonUtil.memoryFootprint()
                    }
                } else if (strict) {
                    throw new FileNotFoundException(fileName)
                } else
                    println "[${new Date()} SampleCollection] WARNING: File $fileName not found, skipping"
            }
        }

        this.metadataTable = metadataTable
    }

    /**
     * Internal util to load sample according to software
     *
     * if store is true, will store the sample to sample map
     * otherwise just returns it
     */
    private Sample loadSample(String sampleId, boolean store) {
        Sample sample

        if (lazy && !loadedSamples.contains(sampleId)) {
            println "[${new Date()} SampleCollection] Loading sample $sampleId"

            sample = new Sample(sampleMap[sampleId].sampleMetadata, new LinkedList<Clonotype>())

            filesBySample[sampleId].each { fileName ->
                sample.clonotypes.addAll(SampleUtil.loadClonotypes(fileName, software))
            }

            println "[${new Date()} SampleCollection] Sample loaded, ${sample.clonotypes.size()} clonotypes. " +
                    CommonUtil.memoryFootprint()

            if (store) {
                loadedSamples.add(sampleId)
                sampleMap.put(sampleId, sample)
            }
        } else
            sample = sampleMap[sampleId]

        sample
    }

    /**
     * Lists all unique sample pairs in a given collection.
     * Pairs (i, j) are chosen such as j > i, no (i, i) pairs allowed.
     * @return a list of sample pairs
     */
    List<SamplePair> listPairs() {
        def samplePairs = new ArrayList()

        // Lazy load all samples
        sampleMap.keySet().each { loadSample(it, true) }

        for (int i = 0; i < size(); i++)
            for (int j = i + 1; j < size(); j++)
                samplePairs.add(this[i, j])

        samplePairs
    }

    /**
     * Get sample by id
     * @param sampleId sample id
     * @return sample
     */
    Sample getAt(String sampleId) {
        sampleMap[sampleId]
    }

    /**
     * Get sample by index, according to current ordering
     * @param i sample index
     * @return sample
     */
    Sample getAt(int i) {
        if (i < 0 || i >= metadataTable.sampleCount)
            throw new IndexOutOfBoundsException()

        sampleMap[metadataTable.getRow(i).sampleId]
    }

    /**
     * Gets sample pair by indices, according to current ordering
     * @param i index of first sample in pair
     * @param j index of second sample in pair
     * @return sample pair
     */
    SamplePair getAt(int i, int j) {
        new SamplePair(this[i], this[j], i, j)
    }

    Iterator iterator() {
        def iter = metadataTable.sampleIterator
        return [
                hasNext: {
                    iter.hasNext()
                },
                next   : {
                    def sampleId = iter.next()
                    loadSample(sampleId, false)
                }] as Iterator
    }

    /**
     * Gets the number of samples
     * @return number of samples
     */
    int size() {
        sampleMap.size()
    }

    @Override
    String toString() {
        "samples=" + this.collect { it.sampleMetadata.sampleId }.join(",") + "\nmetadata=" + metadataTable.join(",")
    }
}
