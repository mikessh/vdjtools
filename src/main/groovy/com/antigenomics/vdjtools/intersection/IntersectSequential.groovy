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

import com.antigenomics.vdjtools.CommonUtil
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleUtil

def cli = new CliBuilder(usage: "IntersectSequential [options] sample1 sample2 sample3 ... output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.c(longOpt: "collapse", argName: "int", args: 1,
        "Generate a collapsed overlap table for visualization purposes with a specified number of top clones.")
cli.m(longOpt: "metadata", argName: "string", args: 1, "Name of tab-delimited metadata file. " +
        "First column should contain sample names, " +
        "second column should contain time course name/units in header and time point values as rows")
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

if (sampleFileNames.unique().size() != sampleFileNames.size())
    println "[ERROR] Duplicate samples not allowed"

def scriptName = getClass().canonicalName.split("\\.")[-1]

def label = "sample"
Map<String, Double> timePointsMap
if (opt.m) {
    timePointsMap = new HashMap<>()
    new File(opt.m.toString()).withReader { reader ->
        label = reader.readLine().split("\t")[1]
        def line
        while ((line = reader.readLine()) != null) {
            def splitLine = line.split("\t")
            timePointsMap.put(sampleFileNames.find { it.contains(splitLine[0]) }, splitLine[1].toDouble())
        }
    }
    if (timePointsMap.size() != sampleFileNames.size()) {
        println "[ERROR] Time points provided (${timePointsMap.values()}) " +
                "don't match the number of samples (n=${sampleFileNames.size()})"
    }
} else {
    timePointsMap = (0..<sampleFileNames.size()).collectEntries { [(sampleFileNames[it]): it] }
}

// Sort samples according to time

timePointsMap = timePointsMap.sort { it.value }

sampleFileNames = timePointsMap.collect { it.key }
def timePoints = timePointsMap.collect { it.value }

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
    // Todo: MDS plot with line or smth else here

    if (opt.c) {
        // println timePoints.join(";")
        CommonUtil.executeR("stack.r", label, timePoints.join(";"),
                outputFilePrefix + "_table_collapsed.txt", outputFilePrefix + "_stackplot.pdf")
    }
}