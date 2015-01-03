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
 * Last modified on 17.11.2014 by mikesh
 */


package com.antigenomics.vdjtools.db

import com.antigenomics.vdjdb.core.db.CdrDatabase
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

// todo: organize CLI for database usage
def cli = new CliBuilder(usage: "ScanDatabase [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.D(longOpt: "database", argName: "string", args: 1, "Path to an external database file.")
cli._(longOpt: "filter", argName: "logical expression(__field__,...)", args: 1,
        "Logical filter on database columns. Supports Regex, .contains(), .startsWith(), etc")
cli.d(longOpt: "details", "Will output detailed DB query files.")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.f(longOpt: "fuzzy", "Will query database allowing at most 2 substitutions, 1 deletion and 1 insertion, " +
        "but no more than 2 mismatches simultaneously.")
cli._(longOpt: "v-match", "V segments are required to match")
cli._(longOpt: "j-match", "V segments are required to match")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(-1)
}

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

def software = Software.byName(opt.S), dbName = (String) (opt.D ?: null), fuzzy = (boolean) opt.f,
    vMatch = (boolean) opt."v-match", jMatch = (boolean) opt."j-match",
    details = (boolean) opt.d, filter = (String) (opt.'filter' ?: null),
    outputFileName = opt.arguments()[-1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) to process"

//
// Annotation
//

def database = dbName ? new CdrDatabase(dbName, filter) : new CdrDatabase(filter)
def databaseBrowser = new DatabaseBrowser(vMatch, jMatch, fuzzy)

println "[${new Date()} $scriptName] Annotating sample(s) & writing results"

new File(formOutputPath(outputFileName, "annot", dbName ?: "default", "summary")).withPrintWriter { pwSummary ->
    def header = "##FILTER=\"$filter\"\n"
    header += "#$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\tdiversity\t" +
            BrowserResult.HEADER

    pwSummary.println(header)

    sampleCollection.eachWithIndex { Sample sample, int ind ->
        def sampleId = sample.sampleMetadata.sampleId
        def browserResult = databaseBrowser.query(sample, database)

        println "[${new Date()} $scriptName] ${ind + 1} sample(s) prepared"

        // Global stats
        pwSummary.println([sampleId, sample.sampleMetadata, sample.diversity, browserResult].join("\t"))

        // Write full summary
        if (details) {
            new File(formOutputPath(outputFileName, "annot", dbName ?: "default", sampleId)).withPrintWriter { pwDetails ->
                pwDetails.println("#" + CdrMatch.HEADER + "\t" + database.ANNOTATION_HEADER)
                browserResult.each { match ->
                    pwDetails.println(match + "\t" + match.subject.annotation.join("\t"))
                }
            }
        }

        println "[${new Date()} $scriptName] ${ind + 1} sample(s) done"
    }
}

println "[${new Date()} $scriptName] Finished"