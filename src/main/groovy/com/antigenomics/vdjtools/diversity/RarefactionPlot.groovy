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
 *
 * Last modified on 18.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.util.RUtil

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath

def STEPS_DEFAULT = "101", I_TYPE_DEFAULT = IntersectionType.Strict
def cli = new CliBuilder(usage: "RarefactionPlot [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")

// general

cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")

// algorithm

cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.s(longOpt: "steps", argName: "int", args: 1, "Number of steps (points) in the rarefaction curve " +
        "(including 0 and the observed diversity). [default=$STEPS_DEFAULT]")

// plotting:

cli.l(longOpt: "label", argName: "string", args: 1,
        "[plotting] Name of metadata column which should be used as label")
cli.f(longOpt: "factor", argName: "string", args: 1,
        "[plotting] Name of metadata column which should be used as a coloring factor")
cli.n(longOpt: "factor-numeric",
        "[plotting] Treat factor values as numeric and use a gradient color scale")
cli._(longOpt: "wide-plot", "[plotting] Will use wide layout for plot")
cli._(longOpt: "label-exact", "[plotting] Will use corresponding sample size for x coordinate of a label. " +
        "Positions all labels at the biggest sample's size if not set.")

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
        println "At least 1 sample file should be provided if not using -m"
    cli.usage()
    System.exit(-1)
}

// Other arguments

def software = Software.byName(opt.S),
    intersectionType = opt.i ? IntersectionType.getByShortName((String) opt.i) : I_TYPE_DEFAULT,
    steps = (opt.s ?: STEPS_DEFAULT).toInteger(),
    optL = opt.'l', optF = opt.'f',
    numericFactor = (boolean) opt.'n',
    widePlot = (boolean) opt.'wide-plot',
    labelExact = (boolean) opt.'label-exact',
    outputPrefix = opt.arguments()[-1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples to analyze"

def sampleStats = sampleCollection.sampleStatistics

//
// Rarefaction analysis
//

def header = ["#$MetadataTable.SAMPLE_ID_COLUMN", sampleCollection.metadataTable.columnHeader,
              Rarefaction.RarefactionPoint.HEADER].flatten().join("\t")

def outputTablePath = formOutputPath(outputPrefix, "rarefaction", intersectionType.shortName)

new File(outputTablePath).withPrintWriter { pw ->
    pw.println(header)

    sampleCollection.eachWithIndex { Sample sample, int i ->
        def sampleId = sample.sampleMetadata.sampleId

        println "[${new Date()} $scriptName] Gathering stats for $sampleId"
        def rarefaction = new Rarefaction(sample, intersectionType)

        println "[${new Date()} $scriptName] Bulding rarefaction curve for $sampleId"
        def rarefactionCurve = rarefaction.build(0, sampleStats.maxCount, steps)

        rarefactionCurve.each {
            pw.println([sampleId, sample.sampleMetadata, it].join("\t"))
        }
    }
}

//
// Plotting for rarefaction analysis
//

println "[${new Date()} $scriptName] Plotting data"

def numeric = RUtil.logical(numericFactor), addLbl = RUtil.logical(optL)

int lblCol = optL ? sampleCollection.metadataTable.getColumnIndex(optL) : 0,
    facCol = optF ? sampleCollection.metadataTable.getColumnIndex(optF) : 0

if (facCol < 0) {
    numeric = RUtil.logical(false)
} else {
    if (numeric == RUtil.logical(true) && sampleCollection.metadataTable.getInfo(optF).numericSamples < 3) {
        println "Switching off numeric option, there were <3 numeric values in corresponding column"
        numeric = RUtil.logical(false)
    }
}

RUtil.execute("rarefaction_curve.r",
        outputTablePath,
        (lblCol + 2).toString(),
        (facCol + 2).toString(),
        numeric,
        addLbl,
        RUtil.logical(widePlot),
        RUtil.logical(labelExact),
        toPlotPath(outputTablePath)
)

println "[${new Date()} $scriptName] Finished"