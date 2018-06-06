package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.graph.DegreeStatistics
import com.antigenomics.vdjtools.graph.DegreeStatisticsCalculator
import com.antigenomics.vdjtools.misc.ExecUtil
import com.antigenomics.vdjtools.misc.StatUtil
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import groovyx.gpars.GParsPool

class DegreeStatisticsAnnotator {
    final DegreeStatisticsCalculator sampleStatistics, controlStatistics

    static final String HEADER = "degree.s\tgroup.count.s\tgroup2.count.s\t" +
            "degree.c\tgroup.count.c\tgroup2.count.c\t" +
            "p.value.g\tp.value.g2"

    DegreeStatisticsAnnotator(DegreeStatisticsCalculator sampleStatistics,
                              DegreeStatisticsCalculator controlStatistics) {
        this.sampleStatistics = sampleStatistics
        this.controlStatistics = controlStatistics
    }

    void annotate(Sample sample) {
        if (sample.annotationHeader) {
            sample.annotationHeader += "\t" + HEADER
            // todo: FIX THREAD COUNT as parameter (?)
            GParsPool.withPool ExecUtil.THREADS, {
                sample.eachParallel { Clonotype clonotype ->
                    def s = sampleStatistics.compute(clonotype),
                        b = controlStatistics.compute(clonotype)
                    clonotype.annotation += "\t" + s + "\t" + b + "\t" +
                            computePValue(s, b) + "\t" + computePValue2(s, b)

                }
            }
        } else {
            sample.annotationHeader = HEADER
            GParsPool.withPool ExecUtil.THREADS, {
                sample.eachParallel { Clonotype clonotype ->
                    def s = sampleStatistics.compute(clonotype),
                        b = controlStatistics.compute(clonotype)
                    clonotype.annotation = s.toString() + "\t" + b + "\t" +
                            computePValue(s, b) + "\t" + computePValue2(s, b)
                }
            }
        }
    }

    static double computePValue(DegreeStatistics sample, DegreeStatistics background) {
        if (sample == DegreeStatistics.UNDEF || background == DegreeStatistics.UNDEF)
            return 1.0

        double p = (background.degree + 1.0d) /
                (background.primaryGroupCount + 1.0d)

        StatUtil.binomialPValue(sample.degree, (double) sample.primaryGroupCount, p) /
                (1.0 - Math.pow(1.0 - p, sample.primaryGroupCount)) // Normalize for conditioning on observing a variant
    }

    static double computePValue2(DegreeStatistics sample, DegreeStatistics background) {
        if (sample == DegreeStatistics.UNDEF || background == DegreeStatistics.UNDEF)
            return 1.0

        double lambda = (sample.secondaryGroupCount + 1.0d) *
                (background.degree + 1.0d) /
                (background.secondaryGroupCount + 1.0d)

        StatUtil.poissonPValue(sample.degree, lambda) /
                (1.0 - Math.exp(-lambda)) // Normalize for conditioning on observing a variant
    }
}
