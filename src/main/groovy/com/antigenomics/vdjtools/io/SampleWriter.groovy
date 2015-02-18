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
 * Last modified on 20.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.io

import com.antigenomics.vdjtools.ClonotypeContainer
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil

import java.util.zip.GZIPOutputStream

/**
 * A class implementing output of Sample and JointSample to plain-text file
 */
public class SampleWriter {
    private final Software software
    private final String header
    private final boolean compress

    /**
     * Creates a sample writer capable to output samples in plain-text files according to specified software format 
     * @param software table layour that will be used during output
     */
    public SampleWriter(Software software) {
        this(software, false)
    }

    /**
     * Creates a sample writer capable to output samples in plain-text files according to specified software format.
     * Will compress the resulting output file if {@code compress = true} and append ".gz" to the file name.
     * @param software table layour that will be used during output
     * @param compress specifies whether to compress resulting output file
     */
    public SampleWriter(Software software, boolean compress) {
        this.software = software
        this.header = (software.headerLineCount > 1 ?
                "$software.name-header-blank\n" * (software.headerLineCount - 1) : "") +
                (software.comment ?: "") +
                software.printFields.join("\t")
        this.compress = compress
    }

    /**
     * Gets buffered writer 
     * @param fileName
     * @return
     */
    private BufferedWriter getWriter(String fileName) {
        if (compress)
            fileName += ".gz"

        def fos = new FileOutputStream(fileName)

        new BufferedWriter(new OutputStreamWriter(compress ?
                new GZIPOutputStream(fos) : fos))
    }

    /**
     * Writes a sample to a given directory/path using conventional sample naming suffix. 
     * Will assume that output to a directory is required if path ends with "/".
     * @param sample sample to write
     * @param outputPrefix output path prefix or output directory (if end with "/")
     */
    public void writeConventional(Sample sample, String outputPrefix) {
        write(sample, ExecUtil.formOutputPath(outputPrefix, sample))
    }

    /**
     * Writes a sample to a given directory/path using conventional sample naming suffix. 
     * Will assume that output to a directory is required if path ends with "/".
     * @param sample sample to write
     * @param outputPrefix output path prefix or output directory (if end with "/")
     * @param top number of top clonotype to write, -1 to write all clonotypes
     * @param collapse specifies whether to store the information on clonotypes that not got it to {@code top}
     *        ones as a separate single entry put at the end of the file
     */
    public void writeConventional(Sample sample, String outputPrefix, int top, boolean collapse) {
        write(sample, ExecUtil.formOutputPath(outputPrefix, sample), top, collapse)
    }

    /**
     * Writes a sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     */
    public void write(ClonotypeContainer sample, String fileName) {
        write(sample, fileName, -1, false)
    }

    /**
     * Writes a sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     * @param top number of top clonotype to write, -1 to write all clonotypes
     * @param collapse specifies whether to store the information on clonotypes that not got it to {@code top}
     *        ones as a separate single entry put at the end of the file
     */
    public void write(ClonotypeContainer sample, String fileName, int top, boolean collapse) {
        def printWriter = getWriter(fileName)

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

        printWriter.close()
    }

    /**
     * Writes a joint sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     * @param top number of top clonotype to write, -1 to write all clonotypes
     * @param collapse specifies whether to store the information on clonotypes that not got it to {@code top}
     *        ones as a separate single entry put at the end of the file
     */
    public void write(JointSample jointSample, String fileName, int top, boolean collapse) {
        def printWriter = getWriter(fileName)
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

        printWriter.close()
    }

    /**
     * Writes a joint sample as a plain-text table to the specified path.
     * @param sample sample to write
     * @param fileName output path
     */
    public void write(JointSample jointSample, String fileName) {
        write(jointSample, fileName, -1, false)
    }
}