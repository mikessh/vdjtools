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
    private final Collection<IntersectMetric> intersectMetrics

    PairedIntersectionSummary(SamplePair samplePair, IntersectionType intersectionType) {
        this(samplePair, intersectionType, IntersectMetric.values())
    }

    PairedIntersectionSummary(SamplePair samplePair,
                              IntersectionType intersectionType,
                              Collection<IntersectMetric> intersectMetrics) {
        this.samplePair = samplePair
        this.jointSample = new JointSample(intersectionType, [samplePair[0], samplePair[1]] as Sample[])
        this.intersectionEvaluator = new IntersectionEvaluator(jointSample)
        this.intersectMetrics = intersectMetrics

        intersectMetrics.each {
            intersectMetricCache.put(it, intersectionEvaluator.computeIntersectionMetric(it))
        }
    }

    public double getDiv1() {
        samplePair[0].diversity
    }

    public int getDiv2() {
        samplePair[1].diversity
    }

    public int getDiv12() {
        jointSample.diversity
    }

    public int getDiv21() {
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
         OUTPUT_FIELDS.collect(), intersectMetrics.collect { it.shortName },
         samplePair[0].sampleMetadata.parent.columnHeader1,
         samplePair[0].sampleMetadata.parent.columnHeader2].flatten().join("\t")
    }

    public String getRow() {
        [samplePair[0].sampleMetadata.sampleId,
         samplePair[1].sampleMetadata.sampleId,
         OUTPUT_FIELDS.collect { this."$it" }, intersectMetrics.collect { intersectMetricCache[it] },
         samplePair[0].sampleMetadata.toString(),
         samplePair[1].sampleMetadata.toString()].flatten().join("\t")
    }
}
