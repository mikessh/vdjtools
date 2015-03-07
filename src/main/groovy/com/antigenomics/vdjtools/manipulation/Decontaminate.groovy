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
 */

package com.antigenomics.vdjtools.manipulation

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.pool.RatioFilter
import com.antigenomics.vdjtools.sample.ClonotypeFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

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
cli._(longOpt: "low-mem", "Will process all sample pairs sequentially, avoiding" +
        " loading all of them into memory. Slower but memory-efficient mode.")
cli.c(longOpt: "compress", "Compress output sample files.")

def opt = cli.parse(args)

if (opt == null) {
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
        println "Output prefix should be provided in case of -m"
    else
        println "At least 2 sample files should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), outputFilePrefix = opt.arguments()[-1],
// todo: implement low-mem version
// this will lead eventually to the problem of removing Clonotype ref from dynamic clonotype
    lowMem = (boolean) opt.'low-mem',
    ratio = (opt.r ?: DEFAULT_CONT_RATIO).toDouble(),
    compress = (boolean) opt.c

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software, false, true) :
        new SampleCollection(opt.arguments()[0..-2], software, false, true)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples prepared"

//
// Vamos a generar un sample aunado
//

println "[${new Date()} $scriptName] Creating sample pool for filtering"

def ratioFilter = new RatioFilter(sampleCollection, ratio)

//
// Go through all sample once more and perform freq-based filtering
//
def sw = new SampleWriter(software, compress)

new File(formOutputPath(outputFilePrefix, "dec", "summary")).withPrintWriter { pw ->
    def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)

    sampleCollection.each { sample ->
        def sampleId = sample.sampleMetadata.sampleId
        println "[${new Date()} $scriptName] Filtering $sampleId"
        def newSample = new Sample(sample, ratioFilter)

        sw.writeConventional(newSample, outputFilePrefix)

        def stats = ratioFilter.getStatsAndFlush()

        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))
    }
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress, "dec:$ratio")

println "[${new Date()} $scriptName] Finished"

