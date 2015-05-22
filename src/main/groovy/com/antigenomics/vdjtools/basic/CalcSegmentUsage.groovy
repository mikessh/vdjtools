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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.util.RUtil

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath

def cli = new CliBuilder(usage: "CalcSegmentUsage [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")
cli._(longOpt: "plot-type", argName: "<pdf|png>", args: 1, "Plot output format [default=pdf]")
cli.p(longOpt: "plot", "Plot V/J usage heatmaps and perform cluster analysis")
cli.n(longOpt: "num-factor", "Numeric factor variable")
cli.l(longOpt: "label", argName: "string", args: 1, "Metadata entry used to annotate the heatmap")
cli.f(longOpt: "factor", argName: "string", args: 1, "Metadata entry used to color samples in the heatmap")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(-1)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

def outputFilePrefix = opt.arguments()[-1],
    unweighted = opt.u,
    plot = (boolean) opt.p,
    plotType = (opt.'plot-type' ?: "pdf").toString()

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])
def metadataTable = sampleCollection.metadataTable

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Compute and output diversity measures, spectratype, etc
//

def segmentUsage = new SegmentUsage(sampleCollection, unweighted)

def outputPathV = formOutputPath(outputFilePrefix, "segments", unweighted ? "unwt" : "wt", "V"),
    outputPathJ = formOutputPath(outputFilePrefix, "segments", unweighted ? "unwt" : "wt", "J")
new File(outputPathV).withPrintWriter { pwV ->
    new File(outputPathJ).withPrintWriter { pwJ ->
        def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" + sampleCollection.metadataTable.columnHeader

        pwV.println(header + "\t" + segmentUsage.vUsageHeader().join("\t"))
        pwJ.println(header + "\t" + segmentUsage.jUsageHeader().join("\t"))

        metadataTable.sampleIterator.each { String sampleId ->
            def sampleString = [sampleId, metadataTable.getRow(sampleId)].join("\t")
            pwV.println(sampleString + "\t" + segmentUsage.vUsageVector(sampleId).collect().join("\t"))
            pwJ.println(sampleString + "\t" + segmentUsage.jUsageVector(sampleId).collect().join("\t"))
        }
    }
}


if (plot) {
    RUtil.execute("vexpr_plot.r",
            outputPathV,
            segmentUsage.vUsageHeader().length.toString(),
            opt.l ? (metadataTable.getColumnIndex(opt.l) + 2).toString() : "0", // first column is sample id
            opt.f ? (metadataTable.getColumnIndex(opt.f) + 2).toString() : "0",
            opt.n ? "TRUE" : "FALSE",
            toPlotPath(outputPathV, plotType)
    )

    RUtil.execute("vexpr_plot.r",
            outputPathJ,
            segmentUsage.jUsageHeader().length.toString(),
            opt.l ? (metadataTable.getColumnIndex(opt.l) + 2).toString() : "0",
            opt.f ? (metadataTable.getColumnIndex(opt.f) + 2).toString() : "0",
            opt.n ? "TRUE" : "FALSE",
            toPlotPath(outputPathJ, plotType)
    )
}

println "[${new Date()} $scriptName] Finished"

