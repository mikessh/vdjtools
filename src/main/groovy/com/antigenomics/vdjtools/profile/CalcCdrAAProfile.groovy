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

package com.antigenomics.vdjtools.profile

import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.util.RUtil

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath


def DEFAULT_AA_PROPERTIES = BasicAminoAcidProperties.INSTANCE.propertyNames.join(","),
    DEFAULT_AA_PROPERTY = "count",
    DEFAULT_BINNING = "V-germ:1,VJ-junc:1,J-germ:1"

def cli = new CliBuilder(usage: "CalcCdrAAProfile [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.u(longOpt: "unweighted", "Will count each clonotype only once. " +
        "[default = clonotype frequency is used for weighting].")
cli.o(longOpt: "property-list", argName: "group1,...", args: 1,
        "Comma-separated list of amino-acid properties to analyze. " +
                "Allowed values: $DEFAULT_AA_PROPERTIES. " +
                "[default = $DEFAULT_AA_PROPERTY]")
cli.r(longOpt: "region-list", argName: "segment1:nbins1,...", args: 1,
        "List of segments to analyze and corresponding bin counts. " +
                "Allowed segments: ${KnownCdr3Regions.INSTANCE.regionNames.join(",")}. " +
                "[default = $DEFAULT_BINNING]")
cli.p(longOpt: "plot", "Plot amino acid property distributions for a specified list of segments.")
cli.f(longOpt: "factor", argName: "string", args: 1, "Metadata entry used to group samples in plot.")
cli._(longOpt: "plot-normalized", "Will normalize regions by the total number of AAs in them.")
cli._(longOpt: "plot-type", argName: "pdf|png", args: 1, "Plot output format [default=pdf]")
cli._(longOpt: "include-cfw", "Consider first and last AAs of CDR3, which are normally conserved C and F/W")


def opt = cli.parse(args)

if (opt == null)
    System.exit(2)

if (opt.h || opt.arguments().size() == 0) {
    cli.usage()
    System.exit(2)
}

// Check if metadata is provided

def metadataFileName = opt.m

if (metadataFileName ? opt.arguments().size() != 1 : opt.arguments().size() < 2) {
    if (metadataFileName)
        println "Only output prefix should be provided in case of -m"
    else
        println "At least 1 sample files should be provided if not using -m"
    cli.usage()
    System.exit(2)
}

// Remaining arguments

def outputFilePrefix = opt.arguments()[-1],
    unweighted = (boolean) opt.u,
    binning = (opt.r ?: DEFAULT_BINNING).split(",").collectEntries {
        def split2 = it.split(":")
        [(KnownCdr3Regions.INSTANCE.getByName(split2[0])): split2[1].toInteger()]
    },
    properties = (opt.o ?: DEFAULT_AA_PROPERTY).split(","),
    plot = (boolean) opt.p,
    plotType = (opt.'plot-type' ?: "pdf").toString(),
    includeCFW = (boolean) opt.'include-cfw'

def badProperties = properties.findAll { !BasicAminoAcidProperties.INSTANCE.propertyNames.contains(it) }

if (badProperties.size() > 0) {
    println "[ERROR] Unknown amino acid properties: ${badProperties.join(",")}. " +
            "Allowed values are $DEFAULT_AA_PROPERTIES"
    System.exit(2)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading sample(s)"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} sample(s) prepared"

//
// Compute and output diversity measures, spectratype, etc
//

def profileBuilder = new Cdr3AAProfileBuilder(binning, !unweighted, !includeCFW, properties)

def outputFileName = formOutputPath(outputFilePrefix, "cdr3aa", "profile", (unweighted ? "unwt" : "wt"))

new File(outputFileName).withPrintWriter { pw ->
    def header = "$MetadataTable.SAMPLE_ID_COLUMN\t" +
            sampleCollection.metadataTable.columnHeader + "\t" +
            "cdr3_segment\tbin\tproperty\tvalue\ttotal"

    pw.println(header)

    def sampleCounter = 0

    sampleCollection.each { Sample sample ->
        def profiles = profileBuilder.create(sample)

        println "[${new Date()} $scriptName] ${++sampleCounter} sample(s) processed"

        profiles.each { profileEntry ->
            def segmentName = profileEntry.key.name
            profileEntry.value.bins.each { bin ->
                bin.summary.each {
                    pw.println([sample.sampleMetadata.sampleId, sample.sampleMetadata,
                                segmentName, bin.index, it.key,
                                it.value, bin.total].join("\t"))
                }
            }
        }
    }
}
if (plot) {
    RUtil.execute("cdr3aa_profile.r",
            outputFileName,
            toPlotPath(outputFileName, plotType),
            opt.f ? (sampleCollection.metadataTable.getColumnIndex(opt.f) + 2).toString() : "0",
            RUtil.logical(opt.'plot-normalized')
    )
}

println "[${new Date()} $scriptName] Finished"