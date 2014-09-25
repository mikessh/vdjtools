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
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.intersection.IntersectionUtil
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.timecourse.DynamicClonotype
import com.antigenomics.vdjtools.timecourse.TimeCourse
import com.antigenomics.vdjtools.util.CommonUtil
import org.apache.commons.io.FilenameUtils

/**
 * Base class to store and handle collections of samples in VDJtools
 * Provides loading samples and their metadata
 * The only limitation is that samples should be processed by same software
 */
class SampleCollection implements Iterable<Sample> {
    private final Map<String, Sample> sampleMap = new HashMap<>()
    private final Map<String, String> filesBySample = new HashMap<>()
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
            reportProgress(it)
        }
    }

    /**
     * Builds a sample collection from a pre-defined list of sample file names.
     * Samples will be assigned to generic metadata table, sample order will be preserved.
     * @param sampleFileNames list of sample file names
     */
    SampleCollection(List<String> sampleFileNames, Software software, boolean lazy) {
        this.software = software
        this.strict = true
        this.lazy = lazy
        this.metadataTable = new MetadataTable()
        sampleFileNames.each { String fileName ->
            def sampleId = FilenameUtils.getBaseName(fileName)
            def sampleMetadata = metadataTable.createRow(sampleId, new ArrayList<String>())
            def clonotypes = lazy ? new ArrayList<Clonotype>() : SampleUtil.loadClonotypes(fileName, software)
            def sample = new Sample(sampleMetadata, clonotypes)
            filesBySample.put(sampleId, fileName)
            sampleMap.put(sampleId, sample)
            reportProgress(sample)
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

        def metadataTable

        new File(sampleMetadataFileName).withReader { reader ->
            metadataTable = new MetadataTable(reader.readLine().split("\t")[2..-1])

            def line, splitLine
            while ((line = reader.readLine()) != null) {
                if (!line.startsWith("#")) {
                    splitLine = line.split("\t")
                    String fileName = splitLine[0], sampleId = splitLine[1]

                    def entries = splitLine.length > 2 ? splitLine[2..-1] : []

                    def inputFile = new File(fileName)

                    if (inputFile.exists() || strict) {
                        def clonotypes = new ArrayList<Clonotype>()

                        if (lazy) {
                            filesBySample.put(sampleId, fileName)
                        } else {
                            clonotypes = SampleUtil.loadClonotypes(fileName, software)
                        }

                        def sampleMetadata = metadataTable.createRow(sampleId, entries)
                        def sample = new Sample(sampleMetadata, clonotypes)
                        sampleMap.put(sampleId, sample)
                        reportProgress(sample)
                    } else
                        println "[${new Date()} SampleCollection] WARNING: File $fileName not found, skipping"
                }
            }
        }

        this.metadataTable = metadataTable
    }

    /**
     * Generates a time course for a given sequential intersection.
     * Only clonotypes met in at least two sequential samples are retained
     * @return clonotype abundance time course
     */
    TimeCourse asTimeCourse() {
        def clonotypeMap = new HashMap<String, Clonotype[]>()

        def intersectionUtil = new IntersectionUtil(IntersectionType.NucleotideV)

        this.eachWithIndex { Sample sample, int sampleIndex ->
            sample.eachWithIndex { Clonotype clonotype, int i ->

                def key = intersectionUtil.generateKey(clonotype)
                def entry = clonotypeMap[key]
                if (!entry)
                    clonotypeMap.put(key, entry = new Clonotype[size()])

                entry[sampleIndex] = clonotype
            }
        }

        new TimeCourse(sampleMap.values() as Sample[], clonotypeMap.values().collect { new DynamicClonotype(it) })
    }

    private int nClonotypes = 0, loadCounter = 0

    private void reportProgress(Sample sample) {
        this.nClonotypes += sample.clonotypes.size()

        if (!lazy) {
            println "[${new Date()} SampleCollection] Loaded ${size()} sample(s) and " +
                    "$nClonotypes clonotypes so far. " + CommonUtil.memoryFootprint()
        } else if (sample.clonotypes.size() > 0) {
            println "[${new Date()} SampleCollection] Loaded sample #${++loadCounter} with " +
                    "${sample.clonotypes.size()} clonotypes so far. " + CommonUtil.memoryFootprint()
        } else {
            println "[${new Date()} SampleCollection] Read ${size()} sample(s)"
        }
    }

    /**
     * Internal util to load sample according to software
     *
     * if store is true, will store the sample to sample map
     * otherwise just returns it
     */
    private Sample lazyLoad(String sampleId, boolean store) {
        Sample sample = sampleMap[sampleId]

        if (lazy && !loadedSamples.contains(sampleId)) {
            def clonotypes = SampleUtil.loadClonotypes(filesBySample[sampleId], software)

            if (store) {
                loadedSamples.add(sampleId)
            } else {
                sample = sample.clone()
                nClonotypes -= clonotypes.size()
            }

            sample.clonotypes.addAll(clonotypes)

            reportProgress(sample)
        }

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
        sampleMap.keySet().each { lazyLoad(it, true) }

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
                    lazyLoad(sampleId, false)
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
