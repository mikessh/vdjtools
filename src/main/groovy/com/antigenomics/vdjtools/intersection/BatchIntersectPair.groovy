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

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.Software

def cli = new CliBuilder(usage: "BatchIntersectPair [options] sample_metadata_file output_prefix")
cli.h("display help message")
cli.I(longOpt: "intersect-type", argName: "string1,string2,..", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$IntersectionType.AminoAcid.shortName' by default.")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

def intersectionTypes

if (opt.I) {
    intersectionTypes = (opt.I as String).split(",").collect {
        def shortName = it.trim()
        def intersectionType = IntersectionType.byName(shortName)
        if (!intersectionType) {
            println "[ERROR] Bad intersection type specified ($shortName). " +
                    "Allowed values are: $IntersectionType.allowedNames"
            System.exit(-1)
        }
    }
} else {
    intersectionTypes = [IntersectionType.AminoAcid]
}

//
// Batch load all samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = new SampleCollection(inputFileName, software, false, false)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Perform intersection for all specified intersection types
//

intersectionTypes.each { IntersectionType intersectionType ->
    println "[${new Date()} $scriptName] Intersecting by $intersectionType"

    def intersectionUtil = new IntersectionUtil(intersectionType)

    def pairedIntersectionMatrix = intersectionUtil.intersectWithinCollection(sampleCollection, false, true)

    println "[${new Date()} $scriptName] Writing results"

    new File(outputFileName + "_" + intersectionType.shortName + ".txt").withPrintWriter { pw ->
        pairedIntersectionMatrix.print(pw)
    }
}