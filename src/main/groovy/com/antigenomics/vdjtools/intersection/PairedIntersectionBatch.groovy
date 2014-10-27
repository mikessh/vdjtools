package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.util.CommonUtil
import com.antigenomics.vdjtools.util.ExecUtil
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger

/**
 * Created by mikesh on 10/26/14.
 */
class PairedIntersectionBatch {
    private final SampleCollection sampleCollection
    private final IntersectionType intersectionType
    private final Collection<IntersectMetric> intersectMetrics
    private final PairedIntersection[][] pairedIntersectionCache
    private final int numberOfSamples

    public PairedIntersectionBatch(SampleCollection sampleCollection,
                                   IntersectionType intersectionType) {
        this(sampleCollection, intersectionType, false, false)
    }

    public PairedIntersectionBatch(SampleCollection sampleCollection,
                                   IntersectionType intersectionType,
                                   boolean store, boolean lowMem) {
        this(sampleCollection, intersectionType, store, lowMem, IntersectMetric.values())
    }

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
                ExecUtil.report(this, "Processed $progr of $totalPairs pairs")
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

            GParsPool.withPool(CommonUtil.THREADS) {
                pairs.eachParallel(intersect)
            }
        }
    }

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

    public String getHeader() {
        ["#sample_id1", "sample_id2",
         PairedIntersection.OUTPUT_FIELDS.collect(), intersectMetrics.collect { it.shortName },
         sampleCollection.metadataTable.columnHeader1,
         sampleCollection.metadataTable.columnHeader2].flatten().join("\t")
    }

    public String getTable() {
        (0..<(numberOfSamples - 1)).collect { int i ->
            ((i + 1)..<numberOfSamples).collect { int j ->
                pairedIntersectionCache[i][j].row
            }
        }.flatten().join("\n")
    }
}