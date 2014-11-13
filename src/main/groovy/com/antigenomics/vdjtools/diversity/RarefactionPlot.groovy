/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 * Last modified on 13.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil
import com.antigenomics.vdjtools.util.RUtil

def STEPS_DEFAULT = "10", RESAMPLES_DEFAULT = "1"
I_TYPE_DEFAULT = IntersectionType.Strict
def cli = new CliBuilder(usage: "RarefactionPlot [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")

cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")

cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")

cli.r(longOpt: "resamples", argName: "int", args: 1,
        "Number of times down-sampling will be performed for each point in rarefaction curve. " +
                "[default=$RESAMPLES_DEFAULT]")
cli.s(longOpt: "steps", argName: "int", args: 1, "Number of steps (points) in the rarefaction curve. " +
        "[default=$RESAMPLES_DEFAULT]")

cli.l(longOpt: "label", argName: "string", args: 1,
        "Name of metadata column which should be used as label")
cli.f(longOpt: "factor", argName: "string", args: 1,
        "Name of metadata column which should be used as a coloring factor")
cli.n(longOpt: "factor-numeric",
        "Treat factor values as numeric and use a gradient color scale")

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
    intersectionType = opt.i ? IntersectionType.byName((String) opt.i) : I_TYPE_DEFAULT,
    steps = (opt.s ?: STEPS_DEFAULT).toInteger(), resamples = (opt.r ?: RESAMPLES_DEFAULT).toInteger(),
    optL = (String) opt.'l', optF = (String) opt.'f', numericFactor = (boolean) opt.'n',
    outputPrefix = opt.arguments()[-1] + ".rarefaction.$intersectionType.shortName"

ExecUtil.ensureDir(outputPrefix)

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples to analyze"

//
// Estimating rarefaction steps
//

def maxCount = 0, minCount = 0

sampleCollection.eachWithIndex { it, ind -> // Sorry we'll have to do it
    maxCount = (int) Math.max(maxCount, it.count)
    minCount = (int) Math.min(minCount, it.count)
}

def rSteps = []

def rStep = minCount / steps

for (int x = 0; x < minCount; x += rStep)
    rSteps.add(x)

rStep = maxCount / steps
for (int x = 0; x < maxCount; x += rStep)
    if (!rSteps.contains(x))
        rSteps.add(x)

//
// Rarefaction analysis
//

def getDiv = { DiversityEstimator diversityEstimator, int x ->
    def subSample = diversityEstimator.downSampler.reSample(x)
    def subSampleDiversityEstimator = new DiversityEstimator(subSample, intersectionType)
    subSampleDiversityEstimator.computeCollapsedSampleDiversity().mean
}

def header = ["#sample_id", sampleCollection.metadataTable.columnHeader,
              "count", "diversity", rSteps].flatten().join("\t")

new File(outputPrefix + ".txt").withPrintWriter { pw ->
    pw.println(header)

    sampleCollection.eachWithIndex { Sample sample, int i ->
        def sampleId = sample.sampleMetadata.sampleId
        def diversityEstimator = new DiversityEstimator(sample, intersectionType)

        def divMax = getDiv(diversityEstimator, (int) sample.count)

        println "[${new Date()} $scriptName] Bulding rarefaction curve for $sampleId"

        for (int k = 0; k < resamples; k++) {
            def y = []

            rSteps.each { int x ->
                if (x > sample.count) {
                    y.add(RUtil.NA)
                } else if (x == 0) {
                    y.add(0)
                } else {
                    y.add(getDiv(diversityEstimator, x))
                }
            }

            pw.println([sampleId, sample.sampleMetadata, sample.count, divMax, y].flatten().join("\t"))
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
    if (sampleCollection.metadataTable.getInfo(optF).numericSamples < 3 && numeric == RUtil.logical(true)) {
        println "Switching off numeric option, there were <3 numeric values in corresponding column"
        numeric = RUtil.logical(false)
    }
}

RUtil.execute("rarefaction_curve.r",
        outputPrefix + ".txt",
        (lblCol + 2).toString(),
        (facCol + 2).toString(),
        numeric,
        addLbl,
        outputPrefix + ".pdf"
)

println "[${new Date()} $scriptName] Finished"