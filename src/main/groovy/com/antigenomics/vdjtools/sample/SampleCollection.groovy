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
import com.antigenomics.vdjtools.ClonotypeUtil
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.Util

class SampleCollection implements Iterable<Sample> {
    private final Map<String, Sample> sampleMap = new HashMap<>()
    private final Map<String, List<String>> filesBySample = new HashMap<>()
    private final HashSet<String> loadedSamples = new HashSet<>()
    private final Software software
    private final boolean strict, lazy
    final List<String> metadataHeader = new ArrayList<>()

    SampleCollection(String sampleMetadataFileName, Software software, boolean strict, boolean lazy) {
        this.software = software
        this.strict = strict
        this.lazy = lazy

        def nSamples = 0, nClonotypes = 0

        new File(sampleMetadataFileName).withReader { reader ->
            metadataHeader.addAll(reader.readLine().split("\t")[2..-1])

            def line, splitLine
            while ((line = reader.readLine()) != null) {
                splitLine = line.split("\t")
                String fileName = splitLine[0], sampleId = splitLine[1]
                def entries = splitLine[2..-1]

                if (entries.size() != metadataHeader.size())
                    throw new Exception("Different number of entries in metadata header and sample $sampleId")

                def inputFile = new File(fileName)

                if (inputFile.exists()) {
                    def clonotypes = new ArrayList<Clonotype>()

                    if (lazy) {
                        def fileList = filesBySample[sampleId]
                        if (!fileList)
                            filesBySample.put(sampleId, fileList = new LinkedList<String>())
                        fileList.add(fileName)
                    } else {
                        clonotypes = ClonotypeUtil.loadClonotypes(fileName, software)
                    }

                    def sample = sampleMap[sampleId]
                    if (!sample) {
                        sampleMap.put(sampleId,
                                new Sample(new SampleMetadata(sampleId, entries),
                                        clonotypes))
                    } else {
                        sample.clonotypes.addAll(clonotypes)
                    }

                    nSamples++

                    if (!lazy) {
                        nClonotypes += clonotypes.size()

                        println "[${new Date()} SampleCollection] Loaded $nSamples samples and " +
                                "$nClonotypes clonotypes so far. " + Util.memoryFootprint()
                    }
                } else if (strict) {
                    throw new FileNotFoundException(fileName)
                } else
                    println "[${new Date()} SampleCollection] WARNING: File $fileName not found, skipping"
            }
        }

        if (metadataHeader[0].startsWith("#"))
            metadataHeader[0] = metadataHeader[0][1..-1]
        metadataHeader.add(0, "#sample_id")
    }

    private Sample loadSample(String sampleId, boolean store) {
        Sample sample

        if (lazy && !loadedSamples.contains(sampleId)) {
            println "[${new Date()} SampleCollection] Loading sample $sampleId"

            sample = new Sample(sampleMap[sampleId].metadata, new LinkedList<Clonotype>())

            filesBySample[sampleId].each { fileName ->
                sample.clonotypes.addAll(ClonotypeUtil.loadClonotypes(fileName, software))
            }

            println "[${new Date()} SampleCollection] Sample loaded, ${sample.clonotypes.size()} clonotypes. " +
                    Util.memoryFootprint()

            if (store) {
                loadedSamples.add(sampleId)
                sampleMap.put(sampleId, sample)
            }
        } else
            sample = sampleMap[sampleId]

        sample
    }

    Collection<SamplePair> listPairs() {
        def samplePairs = new LinkedList()

        // Lazy load all samples
        sampleMap.keySet().each { loadSample(it, true) }

        for (int i = 0; i < sampleMap.values().size(); i++)
            for (int j = i + 1; j < sampleMap.values().size(); j++)
                samplePairs.add(new SamplePair(sampleMap.values()[i], sampleMap.values()[j]))

        samplePairs
    }

    Iterator iterator() {
        def iter = sampleMap.entrySet().iterator()
        return [
                hasNext: {
                    iter.hasNext()
                },
                next   : {
                    def entry = iter.next()
                    loadSample(entry.key, false)
                }] as Iterator
    }

    int size() {
        sampleMap.size()
    }
}
