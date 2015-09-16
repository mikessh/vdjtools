/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */


package com.antigenomics.vdjtools.util

import com.antigenomics.vdjdb.core.db.CdrDatabase
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.util.db.BrowserResult
import com.antigenomics.vdjtools.util.db.CdrMatch
import com.antigenomics.vdjtools.util.db.DatabaseBrowser

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

// todo: organize CLI for database usage
def cli = new CliBuilder(usage: "ScanDatabase [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
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

def dbName = (String) (opt.D ?: null), fuzzy = (boolean) opt.f,
    vMatch = (boolean) opt."v-match", jMatch = (boolean) opt."j-match",
    details = (boolean) opt.d, filter = (String) (opt.'filter' ?: null),
    outputFileName = opt.arguments()[-1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

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