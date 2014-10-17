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

package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil

def cli = new CliBuilder(usage: "CalcBasicStats [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.D(longOpt: "db-name", argName: "string", args: "1",
        "Database name, currently supported: trdb")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.o(longOpt: "one-mismatch", "Will query database allowing a single amino-acid substitution")

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

def software = Software.byName(opt.S), dbName = opt.D ?: "trdb", oneMM = (boolean) opt.o,
    outputFileName = opt.arguments()[-1]

ExecUtil.ensureDir(outputFileName)

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
// Annotation
//

def dbCdrFreqs = new HashMap<String, double[]>()
def database = new CdrDatabase(dbName)

database.each { dbCdrFreqs.put(it, new double[sampleCollection.size()]) }

println "[${new Date()} $scriptName] Annotating sample(s)"

sampleCollection.eachWithIndex { Sample sample, int ind ->
    def sampleAnnotation = new SampleAnnotation(sample, oneMM)
    println "[${new Date()} $scriptName] ${ind + 1} sample(s) prepared"

    sampleAnnotation.getEntryFrequencies(database).each {
        dbCdrFreqs[it.key][ind] += it.value
    }
    println "[${new Date()} $scriptName] ${ind + 1} sample(s) done"
}

//
// Output
//

println "[${new Date()} $scriptName] Writing output"

new File(outputFileName + ".annot.${dbName}.txt").withPrintWriter { pw ->
    pw.println([database.header, sampleCollection.metadataTable.sampleIterator.collect()].flatten().join("\t"))
    dbCdrFreqs.sort { -it.value.collect().sum() }.each {
        def cdr = it.key, values = it.value.collect()
        database.getAnnotationEntries(cdr).each { annotLine ->
            pw.println([cdr, annotLine, values].flatten().join("\t"))
        }
    }
}

println "[${new Date()} $scriptName] Finished"
