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
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection

def N_DEFAULT = "300000", R_STEP_DEFAULT = "250000", R_MAX_DEFAULT = "10000000"
def cli = new CliBuilder(usage: "DiversityStat [options] sample_metadata_file output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.n(longOpt: "num-cells", argName: "integer", args: 1,
        "Number of cells to take for normalized sample diversity estimate. [default = $N_DEFAULT]")
cli.r("Perform rarefaction analysis.")
cli._(longOpt: "r-step", argName: "integer", args: 1,
        "Rarefaction curve step. [default = $R_STEP_DEFAULT]")
cli._(longOpt: "r-max", argName: "integer", args: 1,
        "Rarefaction curve maximum depth. [default = $R_MAX_DEFAULT]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), effectiveSampleSize = (opt.n ?: N_DEFAULT).toInteger(),
    rStep = (opt."r-step" ?: R_STEP_DEFAULT).toInteger(), rMax = (opt."r-max" ?: R_MAX_DEFAULT).toInteger(),
    doRarefaction = opt.r,
    inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

def nResamples = 3, efronDepth = 20, efronCvThreshold = 0.05

def R_EMPTY = ""

//
// Lazy load all samples
//

println "[${new Date()} $scriptName] Reading sample metadata"

def sampleCollection = new SampleCollection(inputFileName, software, false, true)

println "[${new Date()} $scriptName] Processing ${sampleCollection.size()} samples"

//
// Compute and output diversity measures
//

new File(outputFileName + ".diversity.txt").withPrintWriter { pwDiv ->
    def headerDiv = sampleCollection.metadataHeader.join("\t") + "\t" +
            "cells\t" +
            "clones_nt\tclones_aa\t" +
            "cndiv_nt_m\tcndiv_nt_std\t" +
            "efron_nt_m\tefron_nt_std\t" +
            "chao_nt_m\tchao_nt_std\t" +
            "cndiv_aa_m\tcndiv_aa_std\t" +
            "efron_aa_m\tefron_aa_std\t" +
            "chao_aa_m\tchao_aa_std"

    pwDiv.println(headerDiv)

    def rSteps = []
    for (int i = 0; i <= rMax; i += rStep)
        rSteps.add(i)

    def headerR = [sampleCollection.metadataHeader, "cells", "sample_diversity", rSteps].flatten().join("\t")

    new File(outputFileName + ".rarefaction_nt.txt").withPrintWriter { pwRNT ->
        pwRNT.println(headerR)

        new File(outputFileName + ".rarefaction_aa.txt").withPrintWriter { pwRAA ->
            pwRAA.println(headerR)

            int sampleCounter = 0

            sampleCollection.each { Sample sample ->
                def diversityEstimator = new DiversityEstimator(sample)

                // Rarefaction

                if (doRarefaction) {
                    println "[${new Date()} $scriptName] Bulding rarefaction curve"

                    def rNT = [], rAA = []
                    rSteps.each { int i ->
                        if (i > sample.count) {
                            rNT.add(R_EMPTY)
                            rAA.add(R_EMPTY)
                        } else {
                            def ds = diversityEstimator.downSampler.reSample(i)
                            rNT.add(ds.diversityCDR3NT)
                            rAA.add(ds.diversityCDR3AA)
                        }
                    }

                    pwRNT.println([sample.metadata, sample.count, sample.diversityCDR3NT, rNT].flatten().join("\t"))
                    pwRAA.println([sample.metadata, sample.count, sample.diversityCDR3AA, rAA].flatten().join("\t"))
                }

                // Diversity estimates

                println "[${new Date()} $scriptName] Computing diversity estimates"

                def cnDivNT = diversityEstimator.countNormalizedSampleDiversity(effectiveSampleSize, nResamples, false),
                    efronDivNT = diversityEstimator.efronThisted(efronDepth, efronCvThreshold, false),
                    chaoDivNT = diversityEstimator.chao1(false),
                    cnDivAA = diversityEstimator.countNormalizedSampleDiversity(effectiveSampleSize, nResamples, true),
                    efronDivAA = diversityEstimator.efronThisted(efronDepth, efronCvThreshold, true),
                    chaoDivAA = diversityEstimator.chao1(true)

                pwDiv.println([sample.metadata,
                               sample.count,
                               sample.diversityCDR3NT, sample.diversityCDR3AA,
                               cnDivNT, efronDivNT, chaoDivNT,
                               cnDivAA, efronDivAA, chaoDivAA].join("\t"))

                println "[${new Date()} $scriptName] ${++sampleCounter} samples processed"
            }
        }
    }
}