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
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.RUtil

def I_TYPE_DEFAULT = "aa"
def cli = new CliBuilder(usage: "IntersectPair [options] sample1 sample2 output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.c(longOpt: "collapse", argName: "int", args: 1,
        "Generate a collapsed overlap table for visualization purposes with a specified number of top clones.")
cli.p(longOpt: "plot", "Generate a scatterplot to characterize overlapping clonotypes. " +
        "Also generate abundance difference plot if -c option is specified. " +
        "(R installation with ggplot2, grid and gridExtra packages required).")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 3) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), top = (int) ((opt.c ?: "-1").toInteger()),
    sample1FileName = opt.arguments()[0], sample2FileName = opt.arguments()[1],
    outputFilePrefix = opt.arguments()[2]

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Select intersection type

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = IntersectionType.byName(iName)

if (!intersectionType) {
    println "[ERROR] Bad intersection type specified ($iName). " +
            "Allowed values are: $IntersectionType.allowedNames"
    System.exit(-1)
}

//
// Load samples
//

println "[${new Date()} $scriptName] Reading samples $sample1FileName and $sample2FileName"

def sampleCollection = new SampleCollection([sample1FileName, sample2FileName], software, true, false)

//
// Perform an intersection by CDR3NT & V segment
//

println "[${new Date()} $scriptName] Intersecting"

def pairedIntersection = new PairedIntersection(sampleCollection.listPairs()[0], intersectionType, true)
def jointSample = pairedIntersection.jointSample

//
// Generate and write output
//

println "[${new Date()} $scriptName] Writing output"

new File(outputFilePrefix + ".summary.txt").withPrintWriter { pw ->
    pw.println(pairedIntersection.header)
    pw.println(pairedIntersection.row)
}

def sampleWriter = new SampleWriter(software)
sampleWriter.write(jointSample, outputFilePrefix + ".table.txt")

if (top >= 0)
    sampleWriter.write(jointSample, outputFilePrefix + ".table_collapsed.txt", top, true)

if (opt.p) {
    println "[${new Date()} $scriptName] Plotting"

    def sample1 = sampleCollection[0],
        sample2 = sampleCollection[1]

    def log = { double x ->
        Math.log10(x + 1e-7)
    }

    // todo: remake completely
    def xyFile = new File(outputFilePrefix + ".xy.txt")
    xyFile.withPrintWriter { pw ->
        pw.println("x\ty")
        jointSample.each { jointClone ->
            pw.println((0..1).collect { sampleIndex ->
                log(jointClone.getFreq(sampleIndex))
            }.join("\t"))
        }
    }
    xyFile.deleteOnExit()

    def xxFile = new File(outputFilePrefix + ".xx.txt")
    xxFile.withPrintWriter { pw ->
        pw.println("xx")
        sample1.each { pw.println(log(it.freq)) }
    }
    xxFile.deleteOnExit()

    def yyFile = new File(outputFilePrefix + ".yy.txt")
    yyFile.withPrintWriter { pw ->
        pw.println("yy")
        sample2.each { pw.println(log(it.freq)) }
    }
    yyFile.deleteOnExit()

    RUtil.execute("intersect_pair_scatter.r", sample1.sampleMetadata.sampleId, sample2.sampleMetadata.sampleId,
            outputFilePrefix + ".xy.txt", outputFilePrefix + ".xx.txt", outputFilePrefix + ".yy.txt",
            outputFilePrefix + ".scatter.pdf")

    if (opt.c) {
        RUtil.execute("intersect_pair_area.r", sample1.sampleMetadata.sampleId, sample2.sampleMetadata.sampleId,
                outputFilePrefix + ".table_collapsed.txt", outputFilePrefix + ".difference.pdf",
                Math.max(0, software.headerLineCount - 1).toString())
    }
}
