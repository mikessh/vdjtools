package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.graph.DegreeStatisticsCalculator
import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import groovyx.gpars.GParsPool

class DegreeStatisticsAnnotator {
    final DegreeStatisticsCalculator sampleStatistics, controlStatistics

    static final String HEADER = "count.sample\ttotal.sample\tcount.control\ttotal.control"

    DegreeStatisticsAnnotator(DegreeStatisticsCalculator sampleStatistics, DegreeStatisticsCalculator controlStatistics) {
        this.sampleStatistics = sampleStatistics
        this.controlStatistics = controlStatistics
    }

    void annotate(Sample sample) {
        if (sample.annotationHeader) {
            sample.annotationHeader += "\t" + HEADER
            // todo: FIX THREAD COUNT as parameter (?)
            GParsPool.withPool ExecUtil.THREADS, {
                sample.eachParallel { Clonotype clonotype ->
                    clonotype.annotation += "\t" + sampleStatistics.compute(clonotype) + "\t" +
                            controlStatistics.compute(clonotype)
                }
            }
        } else {
            sample.annotationHeader = HEADER
            GParsPool.withPool ExecUtil.THREADS, {
                sample.eachParallel { Clonotype clonotype ->
                    clonotype.annotation = sampleStatistics.compute(clonotype).toString() + "\t" +
                            controlStatistics.compute(clonotype)
                }
            }
        }
    }
}
