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
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleUtil

def cli = new CliBuilder(usage: "IntersectSequence [options] sample1 sample2 sample3 ... output_prefix")
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

if (opt.h || opt.arguments().size() < 4) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S),
    sampleFileNames = opt.arguments()[0..-2],
    outputFilePrefix = opt.arguments()[-1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Load samples
//

println "[${new Date()} $scriptName] Reading samples ${sampleFileNames[0..-2].join(", ")} and ${sampleFileNames[-1]}"

def samples = sampleFileNames.collect { SampleUtil.loadSample(it, software) } as Sample[]

//
// Perform a sequential intersection by CDR3NT & V segment
//

println "[${new Date()} $scriptName] Intersecting"

def sequentialIntersection = new SequentialIntersection(samples, IntersectionType.NucleotideV)

def timeCourse = sequentialIntersection.asTimeCourse()

//
// Generate and write output
//

println "[${new Date()} $scriptName] Writing output"

new File(outputFilePrefix + "_summary.txt").withPrintWriter { pw ->
    // summary statistics: intersection size (count, freq and unique clonotypes)
    // count correlation within intersected set
    pw.println("#" + SequentialIntersection.HEADER)
    pw.println(sequentialIntersection)
}

new File(outputFilePrefix + "_table.txt").withPrintWriter { pw ->
    // all clonotypes in intersection
    timeCourse.print(pw, true)
}

if (opt.c) {
    // top clonotypes in intersection and non-overlapping clonotypes frequency
    int top = (opt.c).toInteger()
    def collapsedTimeCourse = timeCourse.collapseBelow(top)

    new File(outputFilePrefix + "_table_collapsed.txt").withPrintWriter { pw ->
        collapsedTimeCourse.print(pw, true)
    }
}

if (opt.p) {

    if (opt.c) {
        // TODO: stacked area plot
    }
}