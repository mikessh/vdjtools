/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 * Last modified on 7.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.pool.*
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil

def I_TYPE_DEFAULT = "aa"
def cli = new CliBuilder(usage: "PoolSamples [options] " +
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
cli.w(longOpt: "write-cloneset", "Will create a separate file with pooled sample, " +
        "greatly increases memory requirements and potentially creates a very large file.")

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

def software = Software.byName(opt.S), outputPrefix = opt.arguments()[-1],
    writeCloneset = (boolean) opt.'w'

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

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software, false, true) :
        new SampleCollection(opt.arguments()[0..-2], software, false, true)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Pool samples
//

println "[${new Date()} $scriptName] Pooling with $intersectionType .. this may take a while"

def cloneAggrFact = writeCloneset ? new StoringClonotypeAggregatorFactory() : new MaxClonotypeAggregatorFactory()

def sampleAggr = new SampleAggregator(sampleCollection, cloneAggrFact, intersectionType)

//
// Build frequency table for statistics
//

println "[${new Date()} $scriptName] Computing frequency tables for stats"
//def diversityEstimator = new DiversityEstimator(pool)

new File(outputPrefix + ".poolfreqs.txt").withPrintWriter { pwFreq ->
    pwFreq.println("#incidence_count\tread_count")
    sampleAggr.each { MaxClonotypeAggregator cloneAggr ->
        pwFreq.println(cloneAggr.incidenceCount + "\t" + cloneAggr.count)
    }
}

if (writeCloneset) {
    println "[${new Date()} $scriptName] Normalizing and sorting pooled clonotypes. Writing output"
    def writer = new SampleWriter(software)
    writer.write(new PooledSample(sampleAggr), outputPrefix + ".clonetable.txt")
}


