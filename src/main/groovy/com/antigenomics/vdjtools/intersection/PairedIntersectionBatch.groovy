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
 * Last modified on 20.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.util.ExecUtil
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger

/**
 * A class that implements all-vs-all paired intersection for a sample collection and 
 * holds information on those intersections in a manner it could be easily accessed
 */
public class PairedIntersectionBatch {
    private final SampleCollection sampleCollection
    private final IntersectionType intersectionType
    private final Collection<IntersectMetric> intersectMetrics
    private final PairedIntersection[][] pairedIntersectionCache
    private final int numberOfSamples

    /**
     * Intersects clonotype lists for all unique pairs of samples in a given sample collection.
     * Will load both samples into memory for the initialization step. 
     * Pre-computes all intersection metrics.
     * @param sampleCollection a list of samples
     * @param intersectionType clonotype matching rule
     */
    public PairedIntersectionBatch(SampleCollection sampleCollection,
                                   IntersectionType intersectionType) {
        this(sampleCollection, intersectionType, false, false)
    }

    /**
     * Intersects clonotype lists for all unique pairs of samples in a given sample collection.
     * Pre-computes all intersection metrics.
     * @param sampleCollection a list of samples
     * @param intersectionType clonotype matching rule
     * @param store holds all samples in memory if set to {@code true}
     * @param lowMem if set to {@code true}, will not load all samples in memory, but rather load a sample pair at each step
     */
    public PairedIntersectionBatch(SampleCollection sampleCollection,
                                   IntersectionType intersectionType,
                                   boolean store, boolean lowMem) {
        this(sampleCollection, intersectionType, store, lowMem, IntersectMetric.values())
    }

    /**
     * Intersects clonotype lists for all unique pairs of samples in a given sample collection
     * @param sampleCollection a list of samples
     * @param intersectionType clonotype matching rule
     * @param store holds all samples in memory if set to {@code true}
     * @param lowMem if set to {@code true}, will not load all samples in memory, but rather load a sample pair at each step
     * @param intersectMetrics a list of intersection metrics that should be pre-computed
     */
    public PairedIntersectionBatch(SampleCollection sampleCollection,
                                   IntersectionType intersectionType,
                                   boolean store, boolean lowMem,
                                   Collection<IntersectMetric> intersectMetrics) {
        if (store && lowMem)
            throw new Exception("Isn't it illogical to use 'store' and 'lowMem' options simultaneously?")

        this.sampleCollection = sampleCollection
        this.intersectionType = intersectionType
        this.intersectMetrics = intersectMetrics
        this.numberOfSamples = sampleCollection.size()
        this.pairedIntersectionCache = new PairedIntersection[numberOfSamples][numberOfSamples]

        int totalPairs = numberOfSamples * (numberOfSamples - 1) / 2
        def progressCounter = new AtomicInteger()

        def intersect = { SamplePair pair ->
            pairedIntersectionCache[pair.i][pair.j] =
                    new PairedIntersection(pair, intersectionType, store, intersectMetrics)
            int progr
            if ((progr = progressCounter.incrementAndGet()) % 10 == 0) {
                ExecUtil.report(this, "Processed $progr of $totalPairs pairs. " + ExecUtil.memoryFootprint())
            }
        }

        ExecUtil.report(this, "Started batch intersection for $numberOfSamples samples ($totalPairs pairs)")

        if (lowMem) {
            for (int i = 0; i < numberOfSamples - 1; i++) {
                def pairs = sampleCollection.listPairs(i)
                pairs.each(intersect)
            }
        } else {
            def pairs = sampleCollection.listPairs()

            GParsPool.withPool ExecUtil.THREADS, {
                pairs.eachParallel(intersect)
            }
        }
    }

    /**
     * Gets a paired intersection for a given pair of samples 
     * @param i first sample index
     * @param j second sample index
     * @return {@code PairedIntersection} for samples ordered as {@code [i , j]}. Will return {@code null} if {@code i == j}
     * @throws {@code IndexOutOfBoundsException} if {@code i} or {@code j} are not in {@code [0 , numberOfSamples)}
     */
    public PairedIntersection getAt(int i, int j) {
        if (i == j)
            return null // todo: re-implement with dummy batch intersection

        boolean reverse = i > j
        if (reverse)
            (i, j) = [j, i]

        if (i >= numberOfSamples || j < 0)
            throw new IndexOutOfBoundsException()

        reverse ? pairedIntersectionCache[i][j] : pairedIntersectionCache[i][j].reverse
    }

    /**
     * Header string, used for tabular output
     */
    public String getHeader() {
        ["#1_$MetadataTable.SAMPLE_ID_COLUMN", "2_$MetadataTable.SAMPLE_ID_COLUMN",
         PairedIntersection.OUTPUT_FIELDS.collect(), intersectMetrics.collect { it.shortName },
         sampleCollection.metadataTable.columnHeader1,
         sampleCollection.metadataTable.columnHeader2].flatten().join("\t")
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        (0..<(numberOfSamples - 1)).collect { int i ->
            ((i + 1)..<numberOfSamples).collect { int j ->
                pairedIntersectionCache[i][j].toString()
            }
        }.flatten().join("\n")
    }
}