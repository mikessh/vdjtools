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

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.ClonotypeContainer
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil

/**
 * Base class for providing output of Sample and JointSample
 */
class SampleWriter {
    private final Software software
    private final String header

    SampleWriter(Software software) {
        this.software = software
        this.header = (software.headerLineCount > 1 ?
                "$software.name-header-blank\n" * (software.headerLineCount - 1) : "") +
                software.printFields.join("\t")
    }

    public void writeConventional(Sample sample, String outputPrefix) {
        write(sample, ExecUtil.formOutputPath(outputPrefix, sample))
    }

    public void writeConventional(Sample sample, String outputPrefix, int top, boolean collapse) {
        write(sample, ExecUtil.formOutputPath(outputPrefix, sample), top, collapse)
    }

    public void write(ClonotypeContainer sample, String fileName) {
        write(sample, fileName, -1, false)
    }

    public void write(ClonotypeContainer sample, String fileName, int top, boolean collapse) {
        new File(fileName).withPrintWriter { printWriter ->
            top = top > sample.diversity || top < 0 ? sample.diversity : top
            printWriter.println(header)

            long count = 0
            double freq = 0.0

            for (int i = 0; i < top; i++) {
                def clonotype = sample[i]

                if (collapse) {
                    count += clonotype.count
                    freq += clonotype.freq
                }

                printWriter.println(software.printFields.collect {
                    clonotype."$it"
                }.join("\t"))
            }

            if (collapse && top < sample.diversity) {
                // Collapsed
                printWriter.println(software.printFields.collect {
                    if (it == "count")
                        sample.count - count
                    else if (it == "freq")
                        sample.freq - freq
                    else
                        "NotShown"
                }.join("\t"))
            }
        }
    }

    public void write(JointSample jointSample, String fileName, int top, boolean collapse) {
        new File(fileName).withPrintWriter { printWriter ->
            top = top > jointSample.diversity || top < 0 ? jointSample.diversity : top

            def sampleIndices = (0..<jointSample.numberOfSamples)

            printWriter.println(header + "\tpeak\t" +
                    sampleIndices.collect { jointSample.getSample(it).sampleMetadata.sampleId }.join("\t"))

            double collapsedMeanFreq = 0.0
            double[] freqArr = new double[jointSample.numberOfSamples]

            for (int i = 0; i < top; i++) {
                def jointClonotype = jointSample[i],
                    clonotype = jointClonotype.representative

                if (collapse) {
                    collapsedMeanFreq += jointClonotype.baseFreq
                    sampleIndices.each { int j ->
                        freqArr[j] += jointClonotype.getFreq(j)
                    }
                }

                printWriter.println(
                        [software.printFields.collect {
                            if (it == "count")
                                jointClonotype.count
                            else if (it == "freq")
                                jointClonotype.freq
                            else
                                clonotype."$it"
                        },
                         jointClonotype.peak,
                         sampleIndices.collect { int j ->
                             jointClonotype.getFreq(j)
                         }].flatten().join("\t"))
            }

            if (collapse && top < jointSample.diversity) {
                // Collapsed
                def collapsedFreq = jointSample.totalMeanFreq - collapsedMeanFreq,
                    collapsedFreqArr = sampleIndices.collect { int j ->
                        jointSample.getIntersectionFreq(j) - freqArr[j]
                    }
                printWriter.println(
                        [software.printFields.collect {
                            if (it == "count")
                                jointSample.calcCount(collapsedFreq)
                            else if (it == "freq")
                                jointSample.calcFreq(collapsedFreq)
                            else
                                "NotShown"
                        },
                         collapsedFreqArr.findIndexOf { it == collapsedFreqArr.max() },
                         collapsedFreqArr
                        ].flatten().join("\t"))
            }

            // Not in the overlap
            def nonOverlappingFreqArr = sampleIndices.collect { int j ->
                jointSample.getSample(j).freq - jointSample.getIntersectionFreq(j)
            }
            printWriter.println(
                    [software.printFields.collect {
                        if (it == "count")
                            0
                        else if (it == "freq")
                            0.0
                        else
                            "NonOverlapping"
                    },
                     nonOverlappingFreqArr.findIndexOf { it == nonOverlappingFreqArr.max() },
                     nonOverlappingFreqArr
                    ].flatten().join("\t"))
        }
    }

    public void write(JointSample jointSample, String fileName) {
        write(jointSample, fileName, -1, false)
    }
}
