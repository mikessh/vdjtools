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
import com.antigenomics.vdjtools.join.JointSample
import com.antigenomics.vdjtools.join.OccurenceJoinFilter
import com.antigenomics.vdjtools.join.SampleSpecificJoinFilter
import com.antigenomics.vdjtools.parser.SampleWriter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.RUtil

def I_TYPE_DEFAULT = "strict"
def cli = new CliBuilder(usage: "IntersectSequential [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.x(longOpt: "track-sample", argName: "int", args: 1,
        "A zero-based index of time point to track. " +
                "Will otherwise consider all clonotypes that were detected in 2+ samples")
cli.c(longOpt: "collapse", argName: "int", args: 1,
        "Generate a collapsed overlap table for visualization purposes with a specified number of top clones.")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata." +
                "If column named 'time' is present, it will be used to specify time point sequence.")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.t(longOpt: "time", argName: "[t1,t2,t3,...]", args: 1,
        "Time point sequence. Unused if -m is specified.")
cli._(longOpt: "time-label", argName: "string", args: 1,
        "Optional time axis label for plotting.")
// todo: below & at least 3,4,.. samples
cli.p(longOpt: "plot", "Generate a scatterplot to characterize overlapping clonotypes. " +
        "Also generate abundance difference plot if -c option is specified. " +
        "(R installation with ggplot2, grid and gridExtra packages required).")

def opt = cli.parse(args)

if (opt == null) {
    System.exit(-1)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 4) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 3 sample files should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

// Select intersection type

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = IntersectionType.byName(iName)

if (!intersectionType) {
    println "[ERROR] Bad intersection type specified ($iName). " +
            "Allowed values are: $IntersectionType.allowedNames"
    System.exit(-1)
}

def software = Software.byName(opt.S),
    trackSample = (opt.x ?: "-1").toInteger(),
    top = (int) ((opt.c ?: "-1").toInteger()), plot = opt.p,
    timeLabel = opt."time-label" ?: "time",
    outputFilePrefix = opt.arguments()[-1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading in all samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software, true, false) :
        new SampleCollection(opt.arguments()[0..-2], software, true, false)

def metadataTable = sampleCollection.metadataTable

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Sort samples
//

def timePoints = opt.t ? (opt.t.toString())[1..-2].split(",").collect() : // [-10,0.5,1,20]
        (0..<sampleCollection.size()).collect { it.toString() } // value provided with argument or uniform

int n = timePoints.size()

if (!metadataTable.containsColumn("time")) {
    // add corresponding column to metadata
    metadataTable.addColumn("time", timePoints)
}

metadataTable.sort("time", false)

// Take sorted values from metadata
timePoints = metadataTable.getColumn("time").collect { it.asNumeric().toString() }

println "[${new Date()} $scriptName] Time points and samples to be processed:\n" +
        "${timePoints.join(",")}\n${metadataTable.sampleIterator.collect().join(",")}"

//
// Join samples
//

println "[${new Date()} $scriptName] Joining samples by ${intersectionType}"

def jointSample = new JointSample(intersectionType, sampleCollection.collect() as Sample[],
        trackSample > -1 ? new SampleSpecificJoinFilter(trackSample) : new OccurenceJoinFilter(2))

//
// Write output tables
//

println "[${new Date()} $scriptName] Writing tabular output"

def sampleWriter = new SampleWriter(software)
sampleWriter.write(jointSample, outputFilePrefix + ".table.txt")

if (top >= 0)
    sampleWriter.write(jointSample, outputFilePrefix + ".table_collapsed.txt", top, true)

//
// Write summary output
//

println "[${new Date()} $scriptName] Writing output"

new File(outputFilePrefix + ".summary.txt").withPrintWriter { pw ->
    pw.println("#1_sample_id\t2_sample_id\tvalue\tmetric\t" +
            sampleCollection.metadataTable.columnHeader1 + "\t" +
            sampleCollection.metadataTable.columnHeader2)

    double[][] freqTable = new double[n][n],
               divTable = new double[n][n],
               countTable = new double[n][n]

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            freqTable[i][j] = jointSample.getIntersectionFreq(i, j)
            countTable[i][j] = jointSample.getIntersectionCount(i, j)
            divTable[i][j] = jointSample.getIntersectionDiv(i, j)
        }
    }

    def procTable = { double[][] table ->
        def flatTable = table.collect { it.collect() }.flatten().findAll { !Double.isNaN(it) }
        def min = 0, max = 1
        if (flatTable.size() > 0) {
            min = flatTable.min()
            max = flatTable.max()
        } else {
            // should not happen, but who knows
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                table[i][j] = (table[i][j] - min) / (max - min)
            }
        }
    }

    def toRstring = { double x ->
        Double.isFinite(x) ? x : "NA"
    }

    procTable(freqTable)
    procTable(divTable)
    procTable(countTable)

    for (int i = 0; i < jointSample.numberOfSamples; i++) {
        for (int j = 0; j < jointSample.numberOfSamples; j++) {
            def ids = [i, j].collect { jointSample.getSample(it).sampleMetadata.sampleId }.join("\t"),
                metadata = [i, j].collect { jointSample.getSample(it).sampleMetadata.toString() }.join("\t")

            if (i != j) {
                pw.println(ids + "\t" + toRstring(freqTable[i][j]) + "\tfrequency\t" + metadata)
                pw.println(ids + "\t" + toRstring(divTable[i][j]) + "\tdiversity\t" + metadata)
                pw.println(ids + "\t" + toRstring(countTable[i][j]) + "\tcount\t" + metadata)
            } else {
                pw.println(ids + "\tNA\tfrequency\t" + metadata)
                pw.println(ids + "\tNA\tdiversity\t" + metadata)
                pw.println(ids + "\tNA\tcount\t" + metadata)
            }
        }
    }
}

//
// Plotting via R
//

println "[${new Date()} $scriptName] Writing plots"

if (plot) {
    // Plot all the heatmaps
    RUtil.execute("sequential_intersect_similarity_map.r",
            outputFilePrefix + ".summary.txt",
            outputFilePrefix + ".summary.pdf")

    if (top >= 0) {
        // Plot a stack plot of top X clonotype abundances
        RUtil.execute("sequential_intersect_stack.r",
                timeLabel,
                timePoints.join(","),
                outputFilePrefix + ".table_collapsed.txt",
                outputFilePrefix + ".stackplot.pdf")

        // Plot a "heatcourse" plot of top X clonotype abundances
        RUtil.execute("sequential_intersect_heatcourse.r",
                timeLabel,
                timePoints.join(","),
                outputFilePrefix + ".table_collapsed.txt",
                outputFilePrefix + ".heatcourse.pdf")
    }
}