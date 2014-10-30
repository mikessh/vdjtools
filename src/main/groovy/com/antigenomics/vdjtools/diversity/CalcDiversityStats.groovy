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

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.intersection.IntersectionUtil
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil
import com.antigenomics.vdjtools.util.RUtil

def N_DEFAULT = "300000", R_STEP_DEFAULT = "100000", R_MAX_DEFAULT = "1000000",
    I_TYPES_DEFAULT = [IntersectionType.AminoAcid, IntersectionType.Nucleotide]
def cli = new CliBuilder(usage: "CalcDiversityStats [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.n(longOpt: "num-cells", argName: "integer", args: 1,
        "Number of cells to take for normalized sample diversity estimate. [default = $N_DEFAULT]")
cli.r("Perform rarefaction analysis.")
cli.i(longOpt: "intersect-type", argName: "string1,string2,..", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '${I_TYPES_DEFAULT.collect { it.shortName }.join(",")}' by default.")
cli._(longOpt: "r-step", argName: "integer", args: 1,
        "Rarefaction curve step. [default = $R_STEP_DEFAULT]")
cli._(longOpt: "r-max", argName: "integer", args: 1,
        "Rarefaction curve maximum depth. [default = $R_MAX_DEFAULT]")
cli.p(longOpt: "plot", "Plots rarefaction curves. " +
        "(R installation with ggplot2 and reshape packages is reuqired).")

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

// Other arguments

def software = Software.byName(opt.S), effectiveSampleSize = (opt.n ?: N_DEFAULT).toInteger(),
    rStep = (opt."r-step" ?: R_STEP_DEFAULT).toInteger(), rMax = (opt."r-max" ?: R_MAX_DEFAULT).toInteger(),
    doRarefaction = opt.r, plot = opt.p,
    outputPrefix = opt.arguments()[-1]

// Build a list of intersection types to apply

def intersectionTypes

if (opt.i) {
    intersectionTypes = (opt.i as String).split(",").collect {
        def shortName = it.trim()
        def intersectionType = IntersectionType.byName(shortName)
        if (!intersectionType) {
            println "[ERROR] Bad intersection type specified ($shortName). " +
                    "Allowed values are: $IntersectionType.allowedNames"
            System.exit(-1)
        }
        intersectionType
    }
} else {
    intersectionTypes = I_TYPES_DEFAULT
}

ExecUtil.ensureDir(outputPrefix)

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Defaults

def nResamples = 3, efronDepth = 20, efronCvThreshold = 0.05

def R_EMPTY = "NA"

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Compute and output diversity measures
//

// for r plot
def rSteps = []
for (int i = 0; i <= rMax; i += rStep)
    rSteps.add(i)

new File(outputPrefix + ".diversity.txt").withPrintWriter { pwDiv ->
    def headerDiv = "#sample_id\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            "cells\t" +
            intersectionTypes.collect { i ->
                ["clones", "normdiv_m", "normdiv_s", "efron_m", "efron_s", "chao_m", "chao_s"].collect { m ->
                    "${m}_$i.shortName"
                }
            }.flatten().join("\t")

    pwDiv.println(headerDiv)

    def headerR = ["#sample_id", sampleCollection.metadataTable.columnHeader,
                   "cells", "sample_diversity", rSteps].flatten().join("\t")


    intersectionTypes.each { IntersectionType i ->
        new File(outputPrefix + ".rarefaction_${i.shortName}.txt").withPrintWriter { pwR ->
            pwR.println(headerR)
        }
    }

    sampleCollection.each { Sample sample ->
        def diversityRow = [sample.sampleMetadata.sampleId, sample.sampleMetadata,
                            sample.count]

        intersectionTypes.each { IntersectionType intersectionType ->
            def diversityEstimator = new DiversityEstimator(sample, intersectionType)
            def basediv = diversityEstimator.computeCollapsedSampleDiversity(),
                normdiv = diversityEstimator.computeNormalizedSampleDiversity(effectiveSampleSize, nResamples),
                efron = diversityEstimator.computeEfronThisted(efronDepth, efronCvThreshold),
                chao = diversityEstimator.computeChao1()


            diversityRow.addAll([basediv.mean,
                                 normdiv.mean, normdiv.std,
                                 efron.mean, efron.std,
                                 chao.mean, chao.std])

            if (doRarefaction) {
                println "[${new Date()} $scriptName] Bulding rarefaction curve for '${intersectionType.shortName}'"

                for (int k = 0; k < nResamples; k++) {
                    def y = []
                    rSteps.each { int x ->
                        if (x > sample.count) {
                            y.add(R_EMPTY)
                        } else {
                            def subSample = diversityEstimator.downSampler.reSample(x)
                            def subSampleDiversityEstimator = new DiversityEstimator(subSample, intersectionType)
                            y.add(subSampleDiversityEstimator.computeCollapsedSampleDiversity().mean)
                        }
                    }

                    new File(outputPrefix + ".rarefaction_${intersectionType.shortName}.txt").withWriterAppend { pwR ->
                        pwR.println([sample.sampleMetadata.sampleId, sample.sampleMetadata,
                                     sample.count, basediv.mean, y].flatten().join("\t"))
                    }
                }
            }
        }

        pwDiv.println(diversityRow.join("\t"))
    }
}

if (plot) {
    println "[${new Date()} $scriptName] Plotting data"

    // todo: label and color selection, consider not allowing

    def datasets = []

    sampleCollection.metadataTable.sampleIterator.each {
        for (int i =0;i<nResamples;i++)
            datasets.add(it)
    }

    intersectionTypes.each { IntersectionType i ->
        RUtil.execute("rarefaction_curve.r",
                datasets.flatten().join(","), rSteps.join(","),
                outputPrefix + ".rarefaction_${i.shortName}.txt",
                outputPrefix + ".rarefaction_${i.shortName}.pdf"
        )
    }
}

println "[${new Date()} $scriptName] Finished"
