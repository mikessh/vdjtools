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

package com.antigenomics.vdjtools.profile

import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def DEFAULT_AA_GROUPS = BasicAminoAcidProperties.INSTANCE.groupNames.join(","),
    DEFAULT_N_BINS = "10"

def cli = new CliBuilder(usage: "CalcCdrAAProfile [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")
cli.l(longOpt: "property-group-list", argName: "group1,group2,...", args: 1,
        "Comma-separated list of amino-acid property groups to analyze. " +
                "Allowed values: $DEFAULT_AA_GROUPS. " +
                "[default = use all]")
cli.n(longOpt: "bins-count", argName: "integer", args: 1,
        "Number of bins in the profile. [default = $DEFAULT_N_BINS]")


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

def outputFilePrefix = opt.arguments()[-1],
    unweighted = (boolean) opt.u,
    nBins = (opt.n ?: DEFAULT_N_BINS).toInteger(),
    propertyGroups = (opt.l ?: DEFAULT_AA_GROUPS).split(",")

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) prepared"

//
// Compute and output diversity measures, spectratype, etc
//

def profileBuilder = new Cdr3AminoAcidProfileBuilder(nBins, !unweighted, propertyGroups)

new File(formOutputPath(outputFilePrefix, "cdr3aa.profile")).withPrintWriter { pw ->
    def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            "bin\tproperty.group\tproperty\tcount\ttotal"

    pw.println(header)

    def sampleCounter = 0

    sampleCollection.each { Sample sample ->
        def profile = profileBuilder.create(sample)

        println "[${new Date()} $scriptName] ${++sampleCounter} sample(s) processed"

        profile.bins.each { bin ->
            bin.summary.each { groupEntry ->
                def groupName = groupEntry.key
                groupEntry.value.each { propertyEntry ->
                    pw.println([sample.sampleMetadata.sampleId, sample.sampleMetadata,
                                bin.index, groupName, propertyEntry.key,
                                propertyEntry.value, bin.total].join("\t"))
                }
            }
        }
    }
}

println "[${new Date()} $scriptName] Finished"