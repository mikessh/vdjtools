package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SamplePair

/**
 * Created by mikesh on 10/26/14.
 */
class PairedIntersectionSummary {
    private final SamplePair samplePair
    private final JointSample sampleJoin
    private final IntersectionUtil intersectionUtil
    private final IntersectionEvaluator intersectionEvaluator
    private final Map<IntersectMetric, Double> intersectMetricCache = new HashMap<>()

    PairedIntersectionSummary(SamplePair samplePair,
                              IntersectionUtil intersectionUtil) {
        this.samplePair = samplePair
        this.intersectionUtil = intersectionUtil
        this.sampleJoin = new JointSample(intersectionUtil, [samplePair[0], samplePair[1]] as Sample[])
        this.intersectionEvaluator = new IntersectionEvaluator(sampleJoin)
    }

    public double getMetricValue(IntersectMetric intersectMetric) {
        def value = intersectMetricCache[intersectMetric]
        if (!value) {
            intersectMetricCache.put(intersectMetric,
                    value = intersectionEvaluator.computeIntersectionMetric(intersectMetric))
        }
        value
    }

    public double getDiv1() {
        samplePair[0].diversity
    }

    public double getDiv2() {
        samplePair[1].diversity
    }

    public double getDiv12() {
        sampleJoin.diversity
    }

    public double getDiv21() {
        sampleJoin.diversity
    }

    public long getCount1() {
        samplePair[0].count
    }

    public long getCount2() {
        samplePair[1].count
    }

    public long getCount12() {
        sampleJoin.getIntersectionCount(0, 1)
    }

    public long getCount21() {
        sampleJoin.getIntersectionCount(1, 0)
    }

    public double getFreq1() {
        samplePair[0].freq
    }

    public double getFreq2() {
        samplePair[1].freq
    }

    public double getFreq12() {
        sampleJoin.getIntersectionFreq(0, 1)
    }

    public double getFreq21() {
        sampleJoin.getIntersectionFreq(1, 0)
    }
}
