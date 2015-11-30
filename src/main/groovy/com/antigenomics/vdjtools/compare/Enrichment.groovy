/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.compare

import com.antigenomics.vdjtools.group.EnrichmentFilter
import com.antigenomics.vdjtools.group.GroupedSample
import com.antigenomics.vdjtools.group.VJLenScheme
import com.antigenomics.vdjtools.io.SampleFileConnection
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.sample.ClonotypeFilter
import com.antigenomics.vdjtools.sample.IntersectionClonotypeFilter
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def scriptName = getClass().canonicalName.split("\\.")[-1]

def cli = new CliBuilder(usage: "ApplySampleAsFilter [options] " +
        "[sample1 sample2 ... if not -m] control_sample output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.c(longOpt: "compress", "Compress output sample files.")

def opt = cli.parse(args)

if (opt == null) {
    //cli.usage()
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(2)
}

// Check if enough arguments are provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 2 : opt.arguments().size() < 3) {
    if (metadataFileName)
        println "Output prefix and control sample should be provided in case of -m"
    else
        println "At least 1 sample, control sample and output path should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

// IO stuff

def controlFileName = opt.arguments()[-2],
    outputFilePrefix = opt.arguments()[-1],
    compress = (boolean) opt.c

// Parameters

def intersectionType = OverlapType.AminoAcidVJ

//
// Load samples
//

println "[${new Date()} $scriptName] Reading input samples & control sample"

def inputSamples = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-3])

def controlSample = SampleFileConnection.load(controlFileName)

//
// Filter samples
//

println "[${new Date()} $scriptName] Filtering control clonotypes from input samples."

def controlFilter = new IntersectionClonotypeFilter(intersectionType, controlSample, true)

def filteredSamples = SampleCollection.fromSampleList(inputSamples.collect { new Sample(it, controlFilter) })

//
// Compute statistics
//

println "[${new Date()} $scriptName] Pre-processing control sample."

def groupedSample = new GroupedSample(new VJLenScheme())
groupedSample.addAll(controlSample)

def enrichmentFilter = new EnrichmentFilter(false, groupedSample)

println "[${new Date()} $scriptName] Computing enrichment statistics for input clonotypes and writing output."

def sw = new SampleWriter(compress)

new File(formOutputPath(outputFilePrefix, "enrich", "summary")).withPrintWriter { pw ->
    def header = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
            filteredSamples.metadataTable.columnHeader + "\t" +
            ClonotypeFilter.ClonotypeFilterStats.HEADER

    pw.println(header)

    filteredSamples.each { Sample sample ->
        // Filter
        def sampleId = sample.sampleMetadata.sampleId

        println "[${new Date()} $scriptName] Analyzing $sampleId sample."
        def filteredSample = new Sample(sample, enrichmentFilter)

        // print filter stats
        def stats = enrichmentFilter.getStatsAndFlush()
        pw.println([sampleId, sample.sampleMetadata, stats].join("\t"))

        // print output
        sw.writeConventional(filteredSample, outputFilePrefix)
    }
}

filteredSamples.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "enrich:$controlSample.sampleMetadata.sampleId")

println "[${new Date()} $scriptName] Finished"