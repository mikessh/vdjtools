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
 * Last modified on 9.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def cli = new CliBuilder(usage: "CalcSpectratype [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata." +
                "If column named 'time' is present, it will be used to specify time point sequence.")
cli.a(longOpt: "amino-acid", "Will use amino-acid CDR3 sequence lengths instead of nucleotide.")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")

def opt = cli.parse(args)

if (opt == null || opt.arguments().size() == 0)
    System.exit(-1)

if (opt.h) {
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

def software = Software.byName(opt.S),
    outputFilePrefix = opt.arguments()[-1],
    aminoAcid = (boolean) opt.a, unweighted = (boolean) opt.u


def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Compute and output diversity measures, spectratype, etc
//

new File(formOutputPath(outputFilePrefix, "spectratype", (aminoAcid ? "aa" : "nt"), (unweighted ? "unwt" : "wt"))).withPrintWriter { pw ->
    def spectratype = new Spectratype(aminoAcid, unweighted)

    def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" + sampleCollection.metadataTable.columnHeader + "\t" + spectratype.HEADER

    pw.println(header)

    def sampleCounter = 0

    sampleCollection.each { Sample sample ->
        spectratype.addAll(sample)

        println "[${new Date()} $scriptName] ${++sampleCounter} samples processed"

        pw.println([sample.sampleMetadata.sampleId, sample.sampleMetadata, spectratype].join("\t"))

        spectratype.clear()
    }
}

println "[${new Date()} $scriptName] Finished"
