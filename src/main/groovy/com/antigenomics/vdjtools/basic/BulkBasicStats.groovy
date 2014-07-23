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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.system.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection

import java.util.concurrent.atomic.AtomicInteger

def cli = new CliBuilder(usage: "DiversityStat [options] sample_metadata_file output_name")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), inputFileName = opt.arguments()[0],
    outputFileName = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Lazy load all samples
//

println "[${new Date()} $scriptName] Reading sample metadata"

def sampleCollection = new SampleCollection(inputFileName, software, false, true)

println "[${new Date()} $scriptName] Processing ${sampleCollection.size()} samples"

//
// Compute and output diversity measures
//

new File(outputFileName).withPrintWriter { pw ->
    def header = sampleCollection.metadataHeader.join("\t") + "\t" + BasicStats.HEADER

    pw.println(header)

    //def results = new LinkedList<>()

    def sampleCounter = new AtomicInteger()

    //GParsPool.withPool Util.THREADS, {
    //results = sampleCollection.collectParallel { Sample sample ->
    sampleCollection.each { Sample sample ->
        def basicStats = new BasicStats(sample)

        println "[${new Date()} $scriptName] ${sampleCounter.incrementAndGet()} samples processed"

        pw.println([sample.metadata, basicStats].join("\t"))
    }
    //}

    //results.each { pw.println(it) }
}


