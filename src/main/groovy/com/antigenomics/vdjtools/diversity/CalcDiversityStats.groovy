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
 * Last modified on 16.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def N_DEFAULT = "100000", I_TYPES_DEFAULT = [IntersectionType.AminoAcidVJ, IntersectionType.Strict]
def cli = new CliBuilder(usage: "CalcDiversityStats [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.n(longOpt: "num-reads", argName: "integer", args: 1,
        "Number of reads to take for normalized sample diversity estimate. [default = $N_DEFAULT]")
cli.i(longOpt: "intersect-type", argName: "string1,string2,..", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '${I_TYPES_DEFAULT.collect { it.shortName }.join(",")}' by default.")

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

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Defaults

def nResamples = 3, efronDepth = 20, efronCvThreshold = 0.05

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName, software) :
        new SampleCollection(opt.arguments()[0..-2], software)

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples to analyze"

//
// Compute and output diversity measures
//
def maxCells = new HashMap<IntersectionType, Long>()
intersectionTypes.each { maxCells.put(it, 0) }

new File(formOutputPath(outputPrefix, "diversity")).withPrintWriter { pwDiv ->
    def headerDiv = "#sample_id\t" +
            sampleCollection.metadataTable.columnHeader + "\treads\t" +
            intersectionTypes.collect { i ->
                ["clones", "normdiv_m", "normdiv_s", "efron_m", "efron_s", "chao_m", "chao_s",
                 "alpha", "beta", "beta_conf"].collect { m ->
                    "${m}_$i.shortName"
                }
            }.flatten().join("\t")

    pwDiv.println(headerDiv)

    sampleCollection.each { Sample sample ->
        println "[${new Date()} $scriptName] Analyzing $sample.sampleMetadata.sampleId"

        def diversityRow = [sample.sampleMetadata.sampleId, sample.sampleMetadata,
                            sample.count]

        intersectionTypes.each { IntersectionType intersectionType ->
            def diversityEstimator = new DiversityEstimator(sample, intersectionType)
            def basediv = diversityEstimator.computeCollapsedSampleDiversity(),
                normdiv = diversityEstimator.computeNormalizedSampleDiversity(effectiveSampleSize, nResamples),
                efron = diversityEstimator.computeEfronThisted(efronDepth, efronCvThreshold),
                chao = diversityEstimator.computeChao1(),
                freqStat = diversityEstimator.computeFrequencyDistributionStats()


            diversityRow.addAll([basediv.mean,
                                 normdiv.mean, normdiv.std,
                                 efron.mean, efron.std,
                                 chao.mean, chao.std,
                                 freqStat.alpha, freqStat.beta, freqStat.betaConf])

            maxCells.put(intersectionType, Math.max((long) maxCells[intersectionType], sample.count))
        }

        pwDiv.println(diversityRow.join("\t"))
    }
}

println "[${new Date()} $scriptName] Finished"
