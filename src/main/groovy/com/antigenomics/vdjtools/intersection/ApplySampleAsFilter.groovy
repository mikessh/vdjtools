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
 * Last modified on 2.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.IntersectionClonotypeFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil

def scriptName = getClass().canonicalName.split("\\.")[-1]

def cli = new CliBuilder(usage: "ApplySampleAsFilter [options] filter_sample sample1 [sample2 ...] output_prefix")
cli.h("display help message")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection type to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$IntersectionType.Strict.shortName' by default.")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.n(longOpt: "negative", "Will report clonotypes present in filter_sample. The default action is to remove them")

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

if (opt.arguments().size() < 3) {
    println "At least 3 arguments should be provided"
    cli.usage()
    System.exit(-1)
}

// IO stuff

def sampleFileNames = opt.arguments()[0..-2],
    outputPrefix = opt.arguments()[-1]

ExecUtil.ensureDir(outputPrefix)

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

def sampleCollection = new SampleCollection(sampleFileNames.flatten() as List<String>, software)
def clonotypeFilter = new IntersectionClonotypeFilter(intersectionType, sampleCollection[0], negative)

//
// Filter samples
//

println "[${new Date()} $scriptName] Filtering (${negative ? "negative" : "positive"}) and writing output"

def sw = new SampleWriter(software)

new File("${outputPrefix}.filtersummary.txt").withPrintWriter { pw ->
    pw.println("#sample_id\tcells_before\tdiversity_before\tcells_after\tdiversity_after")
    (1..<sampleCollection.size()).each { int index ->
        def sample = sampleCollection[index]

        // Filter
        def filteredSample = new Sample(sample, clonotypeFilter)

        println "[${new Date()} $scriptName] Processed ${index + 1}th sample."

        pw.println([sample.sampleMetadata.sampleId,
                           sample.count, sample.diversity,
                           filteredSample.count, filteredSample.diversity].join("\t"))

        // print output
        sw.write(filteredSample, "${outputPrefix}_${sample.sampleMetadata.sampleId}.txt")
    }
}

println "[${new Date()} $scriptName] Finished"


