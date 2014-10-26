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

package com.antigenomics.vdjtools.parser

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.sample.Sample

/**
 * Base class for providing output of Sample and JointSample
 */
class SampleWriter {
    private final Software software
    private final String header

    SampleWriter(Software software) {
        this.software = software
        this.header = ("$software.name-header-blank\n" * (software.headerLineCount - 1)) + software.printFields.join("\t")
    }

    public void write(Sample sample, PrintWriter printWriter) {
        write(sample, printWriter, -1, false)
    }

    public void write(Sample sample, PrintWriter printWriter, int top, boolean collapse) {
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
            printWriter.println(software.printFields.collect {
                if (it == "count")
                    sample.count - count
                else if (it == "freq")
                    sample.freq - freq
                else
                    "COLLAPSED"
            }.join("\t"))
        }
    }

    public void write(JointSample jointSample, PrintWriter printWriter, int top, boolean collapse) {
        top = top > jointSample.diversity || top < 0 ? jointSample.diversity : top

        def ii = (0..<jointSample.numberOfSamples)

        printWriter.println(header + "\t" +
                ii.collect { jointSample.getSample(it).sampleMetadata.sampleId }.join("\t"))

        long count = 0
        double freq = 0.0
        double[] freqArr = new double[jointSample.numberOfSamples]

        for (int i = 0; i < top; i++) {
            def jointClonotype = jointSample[i],
                clonotype = jointClonotype.representative

            if (collapse) {
                count += jointClonotype.count
                freq += jointClonotype.freq
                ii.each {
                    freqArr[it] += jointClonotype.getFreq(it)
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
                     ii.collect {
                         jointClonotype.getFreq(it)
                     }].flatten().join("\t"))
        }

        if (collapse && top < jointSample.diversity) {
            printWriter.println(
                    [software.printFields.collect {
                        if (it == "count")
                            jointSample.count - count
                        else if (it == "freq")
                            jointSample.freq - freq
                        else
                            "COLLAPSED"
                    },
                     freqArr[ii]
                    ].flatten().join("\t"))
        }
    }

    public void write(JointSample jointSample, PrintWriter printWriter) {
        write(jointSample, printWriter, -1, false)
    }
}
