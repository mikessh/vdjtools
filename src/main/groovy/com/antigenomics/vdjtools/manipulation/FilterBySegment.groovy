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
 * Last modified on 27.2.2015 by mikesh
 */

package com.antigenomics.vdjtools.manipulation

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.*
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def cli = new CliBuilder(usage: "FilterBySegment [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.e(longOpt: "negative", "Will retain only clonotypes that lack specified V/D/J segments")
cli.c(longOpt: "compress", "Compress output sample files.")

cli.v(longOpt: "v-segments", argName: "v1,v2,...", args: 1, "A comma-separated list of Variable segments")
cli.d(longOpt: "d-segments", argName: "d1,d2,...", args: 1, "A comma-separated list of Diversity segments")
cli.j(longOpt: "j-segments", argName: "j1,j2,...", args: 1, "A comma-separated list of Joining segments")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(-1)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

// Remaining arguments

def software = Software.byName(opt.S),
    outputFilePrefix = opt.arguments()[-1],
    compress = (boolean) opt.c,
    negative = (boolean) opt.e

def vFilter = opt.v ? new VFilter(((String) opt.v).split(",") as String[]) : BlankClonotypeFilter.INSTANCE,
    dFilter = opt.d ? new VFilter(((String) opt.d).split(",") as String[]) : BlankClonotypeFilter.INSTANCE,
    jFilter = opt.j ? new VFilter(((String) opt.j).split(",") as String[]) : BlankClonotypeFilter.INSTANCE

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) loaded"

//
// Iterate over samples & filter
//

def writer = new SampleWriter(software, compress)

new File(formOutputPath(outputFilePrefix, "segfilter", "summary")).withPrintWriter { pw ->
    def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)

    sampleCollection.eachWithIndex { sample, ind ->
        def sampleId = sample.sampleMetadata.sampleId
        println "[${new Date()} $scriptName] Filtering $sampleId.."

        def filter = new CompositeClonotypeFilter(negative, vFilter, dFilter, jFilter)
        def filteredSample = new Sample(sample, filter)

        // print output
        writer.writeConventional(filteredSample, outputFilePrefix)

        def stats = filter.getStatsAndFlush()

        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))
    }
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "segfilter:${negative ? "remove" : "keep"}:${opt.v ?: "."}:${opt.d ?: "."}:${opt.j ?: "."}")

println "[${new Date()} $scriptName] Finished"