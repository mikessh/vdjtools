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
 * Last modified on 18.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.basic.SegmentUsage
import com.antigenomics.vdjtools.basic.Spectratype
import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import com.antigenomics.vdjtools.util.MathUtil
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import sun.reflect.generics.reflectiveObjects.NotImplementedException

/**
 * A helper class to compute various intersection metrics for joint intersection
 */
class IntersectionEvaluator {
    public static boolean VERBOSE = true

    private final JointSample jointSample
    private SegmentUsage segmentUsageCache
    private final Spectratype[] spectratypeCache
    private final Map<String, Double> metricCache = new HashMap<>()

    /**
     * Sets up an instance that can compute overlap metrics for a pre-defined set of sample intersections
     * @param jointSample a result of intersection between sample pair(s)
     */
    public IntersectionEvaluator(JointSample jointSample) {
        this.jointSample = jointSample
        this.spectratypeCache = new Spectratype[jointSample.numberOfSamples]
    }

    /**
     * INTERNAL gets spectratype, for spectratype JSD metrics 
     * @param sampleIndex
     * @return
     */
    private Spectratype getSpectratype(int sampleIndex) {
        if (!spectratypeCache[sampleIndex]) {
            spectratypeCache[sampleIndex] = new Spectratype(jointSample.getSample(sampleIndex),
                    jointSample.intersectionType,
                    false)
        }
        spectratypeCache[sampleIndex]
    }

    /**
     * INTERNAL gets segment usage, for V/J/V+J/V*J JSD metrics
     * @param sampleIndex
     * @return
     */
    private SegmentUsage getSegmentUsage() {
        if (!segmentUsageCache) {
            segmentUsageCache = new SegmentUsage((0..<jointSample.numberOfSamples).collect {
                jointSample.getSample(it)
            } as Sample[], false)
        }
        segmentUsageCache
    }

    /**
     * INTERNAL main routine that calculates a specified intersection metric
     * @param metric metric type
     * @param i index of first sample in pair
     * @param j index of second sample in pair
     * @return
     */
    private double _computeIntersectionMetric(IntersectMetric metric,
                                              int i, int j) {
        ExecUtil.report(this, "Computing $metric", VERBOSE)
        switch (metric) {
            case IntersectMetric.Diversity:
                def div1 = jointSample.getSample(i).diversity,
                        div2 = jointSample.getSample(j).diversity,
                        div12 = jointSample.getIntersectionDiv(i, j)
                return div12 / div1 / div2

            case IntersectMetric.Frequency:
                return Math.sqrt(jointSample.getIntersectionFreq(i, j) * jointSample.getIntersectionFreq(j, i))

            case IntersectMetric.Frequency2:
                double F2 = 0;
                jointSample.each {
                    F2 += Math.sqrt(it.getFreq(i) * it.getFreq(j))
                }
                return F2

            case IntersectMetric.Correlation:
                double R = Double.NaN

                int n = jointSample.getIntersectionDiv(i, j)

                if (n > 2) {
                    def x = new double[n],
                        y = new double[n]
                    int k = 0
                    jointSample.each {
                        if (it.present(i) && it.present(j)) {
                            x[k] = it.getFreq(i)
                            y[k] = it.getFreq(j)
                        }
                        k++
                    }

                    R = new PearsonsCorrelation().correlation(x, y)
                }
                return R

            case IntersectMetric.vJSD:
                return MathUtil.JSD(
                        segmentUsage.vUsageVector(0),
                        segmentUsage.vUsageVector(1))

            case IntersectMetric.vjJSD:
                return MathUtil.JSD(
                        [segmentUsage.vUsageVector(0).collect(), segmentUsage.jUsageVector(0).collect()].flatten() as double[],
                        [segmentUsage.vUsageVector(1).collect(), segmentUsage.jUsageVector(1).collect()].flatten() as double[])

            case IntersectMetric.vj2JSD:
                return MathUtil.JSD(
                        segmentUsage.vjUsageMatrix(0).collect().flatten() as double[],
                        segmentUsage.vjUsageMatrix(1).collect().flatten() as double[])

            case IntersectMetric.sJSD:
                return MathUtil.JSD(getSpectratype(i).histogram,
                        getSpectratype(j).histogram)

            default:
                throw new NotImplementedException()
        }
    }

    /**
     * Computes specified intersection metric for a pair of samples
     * @param metric intersection metric type
     * @param i index of first sample in pair
     * @param j index of second sample in pair
     * @return intersection metric value
     */
    public double computeIntersectionMetric(IntersectMetric metric,
                                            int i, int j) {
        def key = [metric.shortName, i, j].join("_")
        def value = metricCache[key]
        if (!value)
            metricCache.put(key, value = _computeIntersectionMetric(metric, i, j))
        value
    }

    /**
     * Computes specified intersection metric for the first pair of samples
     * @param metric intersection metric type
     * @return intersection metric value
     */
    public double computeIntersectionMetric(IntersectMetric intersectMetric) {
        computeIntersectionMetric(intersectMetric, 0, 1)
    }
}
