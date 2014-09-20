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
import com.antigenomics.vdjtools.sample.SampleUtil
import com.antigenomics.vdjtools.util.RUtil

def cli = new CliBuilder(usage: "IntersectPair [options] sample1 sample2 output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
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

def software = Software.byName(opt.S),
    sample1FileName = opt.arguments()[0], sample2FileName = opt.arguments()[1],
    outputFilePrefix = opt.arguments()[2]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Load samples
//

println "[${new Date()} $scriptName] Reading samples $sample1FileName and $sample2FileName"

def sample1 = SampleUtil.loadSample(sample1FileName, software),
    sample2 = SampleUtil.loadSample(sample2FileName, software)

//
// Perform an intersection by CDR3NT & V segment
//

println "[${new Date()} $scriptName] Intersecting"

def intersectionUtil = new IntersectionUtil(IntersectionType.NucleotideV)

def pairedIntersection = intersectionUtil.generatePairedIntersection(sample1, sample2)

def timeCourse = pairedIntersection.asTimeCourse()

//
// Generate and write output
//

println "[${new Date()} $scriptName] Writing output"

new File(outputFilePrefix + ".summary.txt").withPrintWriter { pw ->
    // summary statistics: intersection size (count, freq and unique clonotypes)
    // count correlation within intersected set
    pw.println("#sample1_file\tsample2_file\t" + PairedIntersection.HEADER)
    pw.println(sample1FileName + "\t" + sample2FileName + "\t" + pairedIntersection)
}

new File(outputFilePrefix + ".table.txt").withPrintWriter { pw ->
    // all clonotypes in intersection
    timeCourse.print(pw, true)
}

if (opt.c) {
    // top clonotypes in intersection and non-overlapping clonotypes frequency
    int top = (opt.c).toInteger()
    def collapsedTimeCourse = timeCourse.collapseBelow(top)

    new File(outputFilePrefix + ".table_collapsed.txt").withPrintWriter { pw ->
        collapsedTimeCourse.print(pw, true)
    }
}

def log = { double x ->
    Math.log10(x + 1e-7)
}

if (opt.p) {
    def xyFile = new File(outputFilePrefix + ".xy.txt")
    xyFile.withPrintWriter { pw ->
        pw.println("x\ty")
        timeCourse.each { pw.println(it.frequencies.collect { log(it) }.join("\t")) }
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
                outputFilePrefix + ".table_collapsed.txt", outputFilePrefix + ".difference.pdf")
    }
}
