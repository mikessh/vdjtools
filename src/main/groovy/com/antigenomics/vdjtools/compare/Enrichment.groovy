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

package com.antigenomics.vdjtools.compare

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.group.EnrichmentFilter
import com.antigenomics.vdjtools.group.GroupedSample
import com.antigenomics.vdjtools.group.VJInsScheme
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

def cli = new CliBuilder(usage: "ApplySampleAsFilter [options] " +
        "[sample1 sample2 ... if not -m] control_sample output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
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
        println "Output prefix and control sample should be provided in case of -m"
    else
        println "At least 1 sample, control sample and output path should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

// IO stuff

def controlFileName = opt.arguments()[-2],
    outputFilePrefix = opt.arguments()[-1],
    compress = (boolean) opt.c

// Parameters

def software = Software.byName(opt.S),
    intersectionType = IntersectionType.AminoAcidVJ

//
// Load samples
//

println "[${new Date()} $scriptName] Reading input samples & control sample"

def inputSamples = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-3], software)

def controlSample = SampleFileConnection.load(controlFileName, software)

//
// Filter samples
//

println "[${new Date()} $scriptName] Filtering control clonotypes from input samples."

def controlFilter = new IntersectionClonotypeFilter(intersectionType, controlSample, true)

def filteredSamples = new SampleCollection(inputSamples.collect { new Sample(it, controlFilter) })

//
// Compute statistics
//

println "[${new Date()} $scriptName] Pre-processing control sample."

def groupedSample = new GroupedSample(new VJInsScheme())
groupedSample.addAll(controlSample)

def enrichmentFilter = new EnrichmentFilter(false, groupedSample)

println "[${new Date()} $scriptName] Computing enrichment statistics for input clonotypes and writing output."

def sw = new SampleWriter(software, compress)

new File(formOutputPath(outputFilePrefix, "enrich", "summary")).withPrintWriter { pw ->
    def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" +
            filteredSamples.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)

    filteredSamples.each { Sample sample ->
        // Filter
        def sampleId = sample.sampleMetadata.sampleId

        println "[${new Date()} $scriptName] Analyzing $sampleId sample."
        def filteredSample = new Sample(sample, enrichmentFilter)

        // print filter stats
        def stats = enrichmentFilter.getStatsAndFlush()
        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))

        // print output
        sw.writeConventional(filteredSample, outputFilePrefix)
    }
}

filteredSamples.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "enrich:$controlSample.sampleMetadata.sampleId")

println "[${new Date()} $scriptName] Finished"