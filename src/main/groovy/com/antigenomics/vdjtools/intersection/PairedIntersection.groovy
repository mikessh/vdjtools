package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SamplePair
import com.antigenomics.vdjtools.util.ExecUtil

/**
 * Created by mikesh on 10/26/14.
 */
class PairedIntersection {
    public static boolean VERBOSE = true

    private final SamplePair samplePair
    private final JointSample jointSample
    private final IntersectionEvaluator intersectionEvaluator
    private final Map<IntersectMetric, Double> intersectMetricCache
    private final Collection<IntersectMetric> intersectMetrics
    public final int div1, div2, div12, div21, count1, count2, count12, count21
    public final double freq1, freq2, freq12, freq21
    private final boolean store
    private final String header1, header2,
                         id1, id2, meta1, meta2

    private PairedIntersection(SamplePair samplePair, JointSample jointSample, IntersectionEvaluator intersectionEvaluator,
                               Collection<IntersectMetric> intersectMetrics, Map<IntersectMetric, Double> intersectMetricCache,
                               int div1, int div2, int div12, int div21,
                               int count1, int count2, int count12, int count21,
                               double freq1, double freq2, double freq12, double freq21,
                               String header1, String header2, String id1, String id2, String meta1, String meta2,
                               boolean store) {
        this.samplePair = samplePair
        this.jointSample = jointSample
        this.intersectionEvaluator = intersectionEvaluator
        this.intersectMetrics = intersectMetrics
        this.intersectMetricCache = intersectMetricCache
        this.div1 = div1
        this.div2 = div2
        this.div12 = div12
        this.div21 = div21
        this.count1 = count1
        this.count2 = count2
        this.count12 = count12
        this.count21 = count21
        this.freq1 = freq1
        this.freq2 = freq2
        this.freq12 = freq12
        this.freq21 = freq21
        this.header1 = header1
        this.header2 = header2
        this.id1 = id1
        this.id2 = id2
        this.meta1 = meta1
        this.meta2 = meta2
        this.store = store
    }

    public PairedIntersection(SamplePair samplePair, IntersectionType intersectionType) {
        this(samplePair, intersectionType, false)
    }

    public PairedIntersection(SamplePair samplePair, IntersectionType intersectionType, boolean store) {
        this(samplePair, intersectionType, store, IntersectMetric.values())
    }

    public PairedIntersection(SamplePair samplePair,
                       IntersectionType intersectionType,
                       boolean store,
                       Collection<IntersectMetric> intersectMetrics) {
        this.store = store
        this.samplePair = store ? samplePair : null
        ExecUtil.report(this, "Intersecting samples #${samplePair.i} and ${samplePair.j}", VERBOSE)
        def jointSample = new JointSample(intersectionType, [samplePair[0], samplePair[1]] as Sample[])
        this.jointSample = store ? jointSample : null
        this.intersectionEvaluator = new IntersectionEvaluator(jointSample)
        this.intersectMetrics = intersectMetrics
        this.intersectMetricCache = new HashMap<>()

        intersectMetrics.each {
            intersectMetricCache.put(it, intersectionEvaluator.computeIntersectionMetric(it))
        }

        this.div1 = samplePair[0].diversity
        this.div2 = samplePair[1].diversity
        this.div12 = jointSample.diversity
        this.div21 = div12
        this.count1 = samplePair[0].count
        this.count2 = samplePair[1].count
        this.count12 = jointSample.getIntersectionCount(0, 1)
        this.count21 = jointSample.getIntersectionCount(1, 0)
        this.freq1 = samplePair[0].freq
        this.freq2 = samplePair[1].freq
        this.freq12 = jointSample.getIntersectionFreq(0, 1)
        this.freq21 = jointSample.getIntersectionFreq(1, 0)

        this.header1 = samplePair[0].sampleMetadata.parent.columnHeader1
        this.header2 = samplePair[0].sampleMetadata.parent.columnHeader2
        this.id1 = samplePair[0].sampleMetadata.sampleId
        this.id2 = samplePair[1].sampleMetadata.sampleId
        this.meta1 = samplePair[0].sampleMetadata.toString()
        this.meta2 = samplePair[1].sampleMetadata.toString()
    }

    public double getMetricValue(IntersectMetric intersectMetric) {
        def value = intersectMetricCache[intersectMetric]
        if (!value) {
            if (!store)
                throw new Exception("Cannot provided value for ${intersectMetric.shortName} as " +
                        "\$store=false and the value is not precomputed")
            else
                intersectMetricCache.put(intersectMetric,
                        value = intersectionEvaluator.computeIntersectionMetric(intersectMetric))
        }
        value
    }

    public Sample getSample1() {
        if (!store)
            throw new Exception("Cannot access this property as \$store=false")
        samplePair[0]
    }

    public Sample getSample2() {
        if (!store)
            throw new Exception("Cannot access this property as \$store=false")
        samplePair[1]
    }

    public JointSample getJointSample() {
        if (!store)
            throw new Exception("Cannot access this property as \$store=false")
        jointSample
    }

    public static final String[] OUTPUT_FIELDS = ["div1", "div2", "div12", "div21",
                                                  "count1", "count2", "count12", "count21",
                                                  "freq1", "freq2", "freq12", "freq21"]

    public PairedIntersection getReverse() {
        new PairedIntersection(store ? samplePair.reverse : null, store ? jointSample.reverse : null,
                intersectionEvaluator, intersectMetrics, intersectMetricCache,
                div2, div1, div21, div12,
                count2, count1, count21, count12,
                freq2, freq1, freq21, freq12,
                header1, header2, id2, id1, meta2, meta1,
                store)
    }

    public String getHeader() {
        ["#1_sample_id", "2_sample_id",
         OUTPUT_FIELDS.collect(), intersectMetrics.collect { it.shortName },
         header1, header2].flatten().join("\t")
    }

    public String getRow() {
        [id1, id2,
         OUTPUT_FIELDS.collect { this."$it" }, intersectMetrics.collect { intersectMetricCache[it] },
         meta1, meta2].flatten().join("\t")
    }
}
