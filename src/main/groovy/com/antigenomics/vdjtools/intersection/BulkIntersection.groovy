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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.Util
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.SamplePair
import groovyx.gpars.GParsPool

def cli = new CliBuilder(usage: "BulkIntersection [options] sample_metadata_file output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Bulk load all samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = new SampleCollection(inputFileName, software)

//
// Do intersection in parallel
//

IntersectionType.values().each { IntersectionType intersectionType ->
    println "[${new Date()} $scriptName] Intersecting by $intersectionType"
    def pairs = sampleCollection.listPairs(), results
    GParsPool.withPool Util.THREADS, {
        results = pairs.collectParallel { SamplePair pair ->
            def intersection = new PairedIntersection(pair.sample1.clonotypes, pair.sample2.clonotypes,
                    intersectionType)
            [pair, intersection.intersect()]
        }
    }

    new File(outputFileName + "_" + intersectionType.shortName + ".txt").withPrintWriter { pw ->
        pw.println("#" +
                [sampleCollection.metadataHeader.collect { "1_$it" },
                 sampleCollection.metadataHeader.collect { "2_$it" },
                 IntersectionResult.HEADER].flatten().join("\t"))
        results.each {
            SamplePair samplePair = it[0]
            IntersectionResult intersectionResult = it[1]
            pw.println(samplePair.toString() + "\t" + intersectionResult.toString())
        }
    }
}