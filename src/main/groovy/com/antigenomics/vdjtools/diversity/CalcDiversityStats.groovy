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

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def I_TYPE_DEFAULT = IntersectionType.Strict, RESAMPLES_DEFAULT = "3"
def cli = new CliBuilder(usage: "CalcDiversityStats [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.x(longOpt: "downsample-to", argName: "integer", args: 1,
        "Number of reads to take for down-sampled sample diversity estimate." +
                "[default=number of reads in the smallest sample]")
cli.X(longOpt: "extrapolate-to", argName: "integer", args: 1,
        "Number of reads to take for extrapolated sample diversity estimate." +
                "[default=number of reads in the largest sample]")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule to apply. Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli._(longOpt: "resample-trials", argName: "integer", args: 1,
        "Number of resamples for corresponding estimator. [default = $RESAMPLES_DEFAULT]")

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

def software = Software.byName(opt.S),
    intersectionType = opt.i ? IntersectionType.getByShortName((String) opt.i) : I_TYPE_DEFAULT,
    resampleCount = (opt."resample-trials" ?: RESAMPLES_DEFAULT).toInteger(),
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

//
// Set up downsample/extrapolate read counts
//

def sampleStats = (!opt.x || !opt.X) ? sampleCollection.sampleStatistics : null
def minReads = (opt.x ?: "$sampleStats.minCount").toInteger(),
    maxReads = (opt.X ?: "$sampleStats.maxCount").toInteger()

//
// Compute and output diversity measures
//

def headerBase = "#$MetadataTable.SAMPLE_ID_COLUMN\t" +
        sampleCollection.metadataTable.columnHeader + "\treads\tdiversity"

def exactOutputPath = formOutputPath(outputPrefix, "diversity", intersectionType.shortName, EstimationMethod.Exact.name),
    resampledOutputPath = formOutputPath(outputPrefix, "diversity", intersectionType.shortName, EstimationMethod.Resampled.name)

new File(exactOutputPath).withPrintWriter { pwExact ->
    new File(resampledOutputPath).withPrintWriter { pwResampling ->
        pwExact.println(headerBase + "\textrapolate_reads\t" + ExactEstimator.HEADER)
        pwResampling.println(headerBase + "\tresample_reads\t" + ResamplingEstimator.HEADER)

        sampleCollection.each { Sample sample ->
            println "[${new Date()} $scriptName] Analyzing $sample.sampleMetadata.sampleId"

            def rowBase = [sample.sampleMetadata.sampleId, sample.sampleMetadata,
                           sample.count, sample.diversity].join("\t")

            def exactEstimator = new ExactEstimator(sample, intersectionType, maxReads),
                resamplingEstimator = new ResamplingEstimator(sample, intersectionType, minReads, resampleCount)

            pwExact.println(rowBase + "\t" + maxReads + "\t" + exactEstimator)
            pwResampling.println(rowBase + "\t" + minReads + "\t" + resamplingEstimator)
        }
    }
}

println "[${new Date()} $scriptName] Finished"
