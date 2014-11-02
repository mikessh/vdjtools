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
import com.antigenomics.vdjtools.basic.SegmentUsage
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil

def I_TYPE_DEFAULT = "aa"
def cli = new CliBuilder(usage: "BatchIntersectPair [options] " +
        "[sample1 sample2 sample3 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli._(longOpt: "low-mem", "Will process all sample pairs sequentially, avoiding" +
        " loading all of them into memory. Slower but memory-efficient mode.")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
// todo: provide an option - redundant (i,j) , (j,i) output for the Großtable
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

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 4) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 3 sample files should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), outputPrefix = opt.arguments()[-1],
    lowMem = (boolean) opt.'low-mem'

ExecUtil.ensureDir(outputPrefix)

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Select intersection type

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = IntersectionType.byName(iName)

if (!intersectionType) {
    println "[ERROR] Bad intersection type specified ($iName). " +
            "Allowed values are: $IntersectionType.allowedNames"
    System.exit(-1)
}

//
// Batch load all samples
//

println "[${new Date()} $scriptName] Reading samples"

boolean store, lazy
(store, lazy) = lowMem ? [false, true] : [true, false]

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software, store, lazy) :
        new SampleCollection(opt.arguments()[0..-2], software, store, lazy)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Perform intersection for all specified intersection types
//

println "[${new Date()} $scriptName] Intersecting by $intersectionType"

//def intersectionUtil = new IntersectionUtil(intersectionType)

//def pairedIntersectionMatrix = intersectionUtil.intersectWithinCollection(sampleCollection, false, true)

SegmentUsage.VERBOSE = false
PairedIntersection.VERBOSE = false
IntersectionEvaluator.VERBOSE = false

def pairedIntersectionBatch = new PairedIntersectionBatch(sampleCollection, intersectionType)

println "[${new Date()} $scriptName] Writing results"

new File(outputPrefix + ".batch_intersect_" + intersectionType.shortName + ".txt").withPrintWriter { pw ->
    pw.println(pairedIntersectionBatch.header)
    pw.println(pairedIntersectionBatch.table)
    //pairedIntersectionMatrix.print(pw)
}