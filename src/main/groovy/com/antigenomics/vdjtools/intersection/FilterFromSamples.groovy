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
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SampleUtil

def scriptName = getClass().canonicalName.split("\\.")[-1]

def cli = new CliBuilder(usage: "FilterFromSamples [options] filter_sample sample1 [sample2 ...] output_prefix")
cli.h("display help message")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection type to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$IntersectionType.AminoAcid.shortName' by default.")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.n(longOpt: "negative", "Will report clonotypes present in filter_sample. The default action is to remove them")

def opt = cli.parse(args)

if (opt == null) {
    //cli.usage()
    System.exit(-1)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

// Check if enough arguments are provided

if (opt.arguments().size() < 3) {
    println "At least 3 arguments should be provided"
    cli.usage()
    System.exit(-1)
}

// IO stuff

def filterSampleFileName = opt.arguments()[0],
    sampleFileNames = opt.arguments()[1..-2],
    outputPrefix = opt.arguments()[-1]

// Parameters

def software = Software.byName(opt.S), intersectionType = IntersectionType.byName((opt.i ?: "aa")),
    negative = opt.n

if (!intersectionType) {
    println "[ERROR] Bad intersection type specified ($opt.i). " +
            "Allowed values are: $IntersectionType.allowedNames"
    System.exit(-1)
}

//
// Load samples
//

println "[${new Date()} $scriptName] Reading samples"

def filterSample = SampleUtil.loadSample(filterSampleFileName, software)
def samples = new SampleCollection(sampleFileNames, software, true)

//
// Filter samples
//

println "[${new Date()} $scriptName] Filtering and writing output"

def intersectionUtil = new IntersectionUtil(intersectionType)

samples.eachWithIndex { Sample sample, int index ->
    // Perform the intersection
    def result = intersectionUtil.generatePairedIntersection(filterSample, sample, true, false)

    // Filter
    Sample newSample

    if (negative) {
        // add only clonotypes in intersection
        newSample = new Sample("filtered_" + sample.sampleMetadata.sampleId, result.clonotypes21)
    } else {
        // remove clonotypes in intersection
        def filterSet = new HashSet<>(result.clonotypes21)
        sample.clonotypes.removeAll { filterSet.contains(it) }
        newSample = sample
    }

    // recompute frequencies and sort
    newSample.renormalize()
    newSample.sort()

    println "[${new Date()} $scriptName] Processed ${index + 1}th sample. " +
            "${negative ? "Retained" : "Removed"} ${(int)(result.freq21 * 10000) / (int)100}% cells " +
            "and $result.div12 clonotypes"

    // print output
    new File(outputPrefix + "_" + sample.sampleMetadata.sampleId + ".txt").withPrintWriter { pw ->
        newSample.print(pw, software, true)
    }
}

println "[${new Date()} $scriptName] Finished"


