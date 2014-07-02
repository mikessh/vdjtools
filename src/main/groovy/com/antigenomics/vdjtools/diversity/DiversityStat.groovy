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
import com.antigenomics.vdjtools.Util
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger


def N_DEFAULT = "300000"
def cli = new CliBuilder(usage: "BulkIntersection [options] sample_metadata_file output_name")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.n(longOpt: "num-cells", argName: "integer", args: 1,
        "Number of cells to take for normalized sample diversity estimate. [default = $N_DEFAULT]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), effectiveSampleSize = (opt.n ?: N_DEFAULT).toInteger(),
    inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

def nResamples = 3, efronDepth = 20, efronCvThreshold = 0.1

//
// Lazy load all samples
//

println "[${new Date()} $scriptName] Reading sample metadata"

def sampleCollection = new SampleCollection(inputFileName, software, false, true)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Compute and output diversity measures
//

new File(outputFileName).withPrintWriter { pw ->
    def header = sampleCollection.metadataHeader.join("\t") + "\t" +
            "cells\tclones_nt\tclones_aa" +
            "cndiv_nt_m\tcndiv_nt_std\t" +
            "efron_nt_m\tefron_nt_std\t" +
            "chao_nt_m\tchao_nt_std" +
            "cndiv_aa_m\tcndiv_aa_std\t" +
            "efron_aa_m\tefron_aa_std\t" +
            "chao_aa_m\tchao_aa_std"

    pw.println(header)

    def results = new LinkedList<>()

    def sampleCounter = new AtomicInteger()

    GParsPool.withPool Util.THREADS, {
        results = sampleCollection.collectParallel { Sample sample ->
            def diversityEstimator = new DiversityEstimator(sample)
            def cnDivNT = diversityEstimator.countNormalizedSampleDiversity(effectiveSampleSize, nResamples, false),
                efronDivNT = diversityEstimator.efronThisted(efronDepth, efronCvThreshold, false),
                chaoDivNT = diversityEstimator.chao1(false),
                cnDivAA = diversityEstimator.countNormalizedSampleDiversity(effectiveSampleSize, nResamples, true),
                efronDivAA = diversityEstimator.efronThisted(efronDepth, efronCvThreshold, true),
                chaoDivAA = diversityEstimator.chao1(true)

            println "[${new Date()} $scriptName] ${sampleCounter.incrementAndGet()} samples processed"

            [sample.metadata, sample.cells, sample.clonotypes,
             cnDivNT, efronDivNT, chaoDivNT,
             cnDivAA, efronDivAA, chaoDivAA].join("\t")
        }
    }

    results.each { pw.println(it) }
}

