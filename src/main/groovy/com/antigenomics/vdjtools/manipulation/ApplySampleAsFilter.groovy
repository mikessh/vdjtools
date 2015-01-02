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
 * Last modified on 2.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.manipulation

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.io.SampleFileConnection
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.ClonotypeFilter
import com.antigenomics.vdjtools.sample.IntersectionClonotypeFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def scriptName = getClass().canonicalName.split("\\.")[-1]

def cli = new CliBuilder(usage: "ApplySampleAsFilter [options] filter_sample " +
        "[sample1 sample2 ... if not -m] filter_sample output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection type to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$IntersectionType.Strict.shortName' by default.")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.n(longOpt: "negative", "Will report clonotypes present in filter_sample. The default action is to remove them")
cli.c(longOpt: "compress", "Compress output sample files.")

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

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 2 : opt.arguments().size() < 3) {
    if (metadataFileName)
        println "Output prefix and filter sample should be provided in case of -m"
    else
        println "At least 1 sample, filter sample and output path should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

// IO stuff

def filterFileName = opt.arguments()[-2],
    compress = (boolean) opt.c,
    outputFilePrefix = opt.arguments()[-1]

// Parameters

def software = Software.byName(opt.S), intersectionType = IntersectionType.byName((opt.i ?: "strict")),
    negative = (boolean) opt.n

if (!intersectionType) {
    println "[ERROR] Bad intersection type specified ($opt.i). " +
            "Allowed values are: $IntersectionType.allowedNames"
    System.exit(-1)
}

//
// Load samples
//

println "[${new Date()} $scriptName] Reading input samples & filter sample"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-3], software)

println "[${new Date()} $scriptName] Loading filter sample"

def filterSample = SampleFileConnection.load(filterFileName, software)

def clonotypeFilter = new IntersectionClonotypeFilter(intersectionType, filterSample, negative)

//
// Filter samples
//

println "[${new Date()} $scriptName] Filtering (${negative ? "negative" : "positive"}) and writing output"

def sw = new SampleWriter(software, compress)

new File(formOutputPath(outputFilePrefix, "asaf", "summary")).withPrintWriter { pw ->
    def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)
    
    sampleCollection.each { Sample sample ->
        // Filter
        def sampleId = sample.sampleMetadata.sampleId
        
        println "[${new Date()} $scriptName] Filtering $sampleId sample."
        def filteredSample = new Sample(sample, clonotypeFilter)

        // print filter stats
        def stats = clonotypeFilter.getStatsAndFlush()
        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))

        // print output
        sw.writeConventional(filteredSample, outputFilePrefix)
    }
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "asaf:$filterSample.sampleMetadata.sampleId:${negative ? "-" : "+"}:$intersectionType.shortName")

println "[${new Date()} $scriptName] Finished"