package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SamplePair

/**
 * Created by mikesh on 10/26/14.
 */
class PairedIntersectionSummary {
    private final SamplePair samplePair
    private final JointSample jointSample
    private final IntersectionEvaluator intersectionEvaluator
    private final Map<IntersectMetric, Double> intersectMetricCache = new HashMap<>()

    PairedIntersectionSummary(SamplePair samplePair,
                              IntersectionType intersectionType) {
        this.samplePair = samplePair
        this.jointSample = new JointSample(intersectionType, [samplePair[0], samplePair[1]] as Sample[])
        this.intersectionEvaluator = new IntersectionEvaluator(jointSample)
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
        jointSample.diversity
    }

    public double getDiv21() {
        jointSample.diversity
    }

    public long getCount1() {
        samplePair[0].count
    }

    public long getCount2() {
        samplePair[1].count
    }

    public long getCount12() {
        jointSample.getIntersectionCount(0, 1)
    }

    public long getCount21() {
        jointSample.getIntersectionCount(1, 0)
    }

    public double getFreq1() {
        samplePair[0].freq
    }

    public double getFreq2() {
        samplePair[1].freq
    }

    public double getFreq12() {
        jointSample.getIntersectionFreq(0, 1)
    }

    public double getFreq21() {
        jointSample.getIntersectionFreq(1, 0)
    }

    JointSample getJointSample() {
        return jointSample
    }

    private static final String[] OUTPUT_FIELDS = ["div1", "div2", "div12", "div21",
                                            "count1", "count2", "count12", "count21",
                                            "freq1", "freq2", "freq12", "freq21"]

    public String getHeader() {
        ["#sample_id1", "sample_id2",
         OUTPUT_FIELDS.collect(), IntersectMetric.values().collect { it.shortName },
         samplePair[0].sampleMetadata.parent.columnHeader1,
         samplePair[0].sampleMetadata.parent.columnHeader2].flatten().join("\t")
    }

    public String getRow() {
        [samplePair[0].sampleMetadata.sampleId,
         samplePair[1].sampleMetadata.sampleId,
         OUTPUT_FIELDS.collect { this."$it" }, IntersectMetric.values().collect { getMetricValue(it) },
         samplePair[0].sampleMetadata.toString(),
         samplePair[1].sampleMetadata.toString()].flatten().join("\t")
    }
}
