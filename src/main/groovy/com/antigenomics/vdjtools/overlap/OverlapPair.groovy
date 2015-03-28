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

package com.antigenomics.vdjtools.overlap

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.RUtil

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath

def I_TYPE_DEFAULT = "strict", TOP_DEFAULT = "20", TOP_MAX = 100
def cli = new CliBuilder(usage: "OverlapPair [options] sample1 sample2 output_prefix")
cli.h("display help message")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.t(longOpt: "top", args: 1, "Number of top clonotypes which will be provided in the collapsed joint table " +
        "and shown on the summary stacked area plot. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot unreadable. [default = $TOP_DEFAULT]")
cli.p(longOpt: "plot", "Generate a scatterplot to characterize overlapping clonotypes. " +
        "Also generate abundance difference plot if -c option is specified. " +
        "(R installation with ggplot2, grid and gridExtra packages required).")
cli.c(longOpt: "compress", "Compress output sample files.")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 3) {
    cli.usage()
    System.exit(-1)
}

def sample1FileName = opt.arguments()[0], sample2FileName = opt.arguments()[1],
    outputPrefix = opt.arguments()[2],
    compress = (boolean) opt.c

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Select overlap type

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = OverlapType.getByShortName(iName)

if (!intersectionType) {
    println "[ERROR] Bad overlap type specified ($iName). " +
            "Allowed values are: $OverlapType.allowedNames"
    System.exit(-1)
}

// Define number of clonotypes to show explicitly

def top = (opt.t ?: TOP_DEFAULT).toInteger()

if (top > TOP_MAX) {
    println "[ERROR] Specified number of top clonotypes should not exceed $TOP_MAX"
    System.exit(-1)
}

//
// Load samples
//

println "[${new Date()} $scriptName] Reading samples $sample1FileName and $sample2FileName"

def sampleCollection = new SampleCollection([sample1FileName, sample2FileName], Software.VDJtools, true, false)

//
// Perform an overlap by CDR3NT & V segment
//

println "[${new Date()} $scriptName] Intersecting"

def pairedIntersection = new Overlap(sampleCollection.listPairs()[0], intersectionType, true)
def jointSample = pairedIntersection.jointSample

//
// Generate and write output
//

println "[${new Date()} $scriptName] Writing output"

new File(formOutputPath(outputPrefix, "paired", intersectionType.shortName, "summary")).withPrintWriter { pw ->
    pw.println(pairedIntersection.header)
    pw.println(pairedIntersection.toString())
}


def sampleWriter = new SampleWriter(compress)
sampleWriter.write(jointSample, formOutputPath(outputPrefix, "paired", intersectionType.shortName, "table"))

def tableCollapsedOutputPath = formOutputPath(outputPrefix, "paired", intersectionType.shortName, "table", "collapsed")
if (top >= 0)
    sampleWriter.write(jointSample, tableCollapsedOutputPath, top, true)

if (opt.p) {
    println "[${new Date()} $scriptName] Plotting"

    def sample1 = sampleCollection[0],
        sample2 = sampleCollection[1]

    def xyFile = new File(outputPrefix + ".xy.txt")
    xyFile.withPrintWriter { pw ->
        pw.println("x\ty")
        jointSample.each { jointClone ->
            pw.println((0..1).collect { sampleIndex ->
                jointClone.getFreq(sampleIndex)
            }.join("\t"))
        }
    }
    xyFile.deleteOnExit()

    def xxFile = new File(outputPrefix + ".xx.txt")
    xxFile.withPrintWriter { pw ->
        pw.println("xx")
        sample1.each { pw.println(it.freq) }
    }
    xxFile.deleteOnExit()

    def yyFile = new File(outputPrefix + ".yy.txt")
    yyFile.withPrintWriter { pw ->
        pw.println("yy")
        sample2.each { pw.println(it.freq) }
    }
    yyFile.deleteOnExit()

    RUtil.execute("intersect_pair_scatter.r", sample1.sampleMetadata.sampleId, sample2.sampleMetadata.sampleId,
            outputPrefix + ".xy.txt", outputPrefix + ".xx.txt", outputPrefix + ".yy.txt",
            formOutputPath(outputPrefix, intersectionType.shortName, "paired", "scatter", "pdf"))

    RUtil.execute("intersect_pair_area.r", sample1.sampleMetadata.sampleId, sample2.sampleMetadata.sampleId,
            tableCollapsedOutputPath, toPlotPath(tableCollapsedOutputPath),
            Math.max(0, Software.VDJtools.headerLineCount - 1).toString())
}
