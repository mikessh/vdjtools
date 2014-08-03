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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.timecourse.DynamicClonotype

def DEFAULT_CONT_RATIO = "20"

def cli = new CliBuilder(usage: "Decontaminate [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.r(longOpt: "ratio", argName: "double", args: 1,
        "Parent-to-child clonotype frequency ratio for contamination filtering [default = $DEFAULT_CONT_RATIO]")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")


def opt = cli.parse(args)

if (opt == null) {
    //cli.usage()
    System.exit(-1)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 3) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 2 sample files should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), outputFilePrefix = opt.arguments()[-1],
    ratio = (opt.r ?: DEFAULT_CONT_RATIO).toDouble()

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software, false, false) :
        new SampleCollection(opt.arguments()[0..-2], software, false)

def nSamples = sampleCollection.size()

def metadataTable = sampleCollection.metadataTable

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Collect dynamic clonotypes
//

println "[${new Date()} $scriptName] Joining samples"

def timeCourse = sampleCollection.asTimeCourse()
timeCourse.sort()

//
// Filter and output
//

def printWriters = metadataTable.sampleIterator.collect { String sampleId ->
    new PrintWriter(new File(outputFilePrefix + "_" + sampleId + ".txt"))
}

printWriters.each { PrintWriter pw ->
    pw.println(software.header)
}

def filteredClonotypes = new int[nSamples],
    filteredFreq = new double[nSamples]

println "[${new Date()} $scriptName] Filtering and writing output"

def clonotypeCounter = 0

timeCourse.each { DynamicClonotype clonotype ->
    def maxFreq = clonotype.frequencies[clonotype.peak]
    (0..<nSamples).each { int sampleIndex ->
        def freq = clonotype.frequencies[sampleIndex]
        if (freq > 0) {
            if (maxFreq < freq * ratio) {
                clonotype[sampleIndex].print(printWriters[sampleIndex], software)
            } else {
                filteredClonotypes[sampleIndex]++
                filteredFreq[sampleIndex] += freq
            }
        }
    }

    if (++clonotypeCounter % 100000 == 0)
        println "[${new Date()} $scriptName] $clonotypeCounter clonotypes processed"
}

new File(outputFilePrefix + "_summary.txt").withPrintWriter { pw ->
    pw.println("#\t" + metadataTable.sampleIterator.collect { "filtered_$it" }.join("\t") + "\t" +
            metadataTable.sampleIterator.collect { "total_$it" }.join("\t"))
    pw.println("frequncy\t" + filteredFreq.collect().join("\t") + "\t" +
            sampleCollection.collect { it.freq }.join("\t"))
    pw.println("clonotypes\t" + filteredClonotypes.collect().join("\t") + "\t" +
            sampleCollection.collect { it.div }.join("\t"))
}

printWriters.each { pw -> pw.close() }

println "[${new Date()} $scriptName] Finished"

