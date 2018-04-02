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

package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.graph.DegreeStatisticsCalculator
import com.antigenomics.vdjtools.graph.DummyClonotypeGroupingFactory
import com.antigenomics.vdjtools.graph.VJClonotypeGroupingFactory
import com.antigenomics.vdjtools.graph.VJLClonotypeGroupingFactory
import com.antigenomics.vdjtools.io.SampleFileConnection
import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.overlap.OverlapType
import com.antigenomics.vdjtools.pool.PooledSample
import com.antigenomics.vdjtools.pool.SampleAggregator
import com.antigenomics.vdjtools.pool.StoringClonotypeAggregatorFactory
import com.antigenomics.vdjtools.sample.SampleCollection


def DEFAULT_SEARCH_SCOPE = "1,0,1"

def cli = new CliBuilder(usage: "CalcDegreeStats [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.b(longOpt: "background", argName: "filename", args: 1,
        "Background (control) sample file to compute random graph statistics. In case not provided," +
                "will pool input samples and use them as control.")
cli.o(longOpt: "search-scope", argName: "s,id,t", args: 1,
        "Search scope: number of substitutions (s), indels (id) and total number of mismatches (t) allowed. " +
                "Default is $DEFAULT_SEARCH_SCOPE")
cli.g(longOpt: "grouping", argName: "type", args: 1,
        "Grouping type: 'dummy', 'vj', 'vjl'. By default will select 'vjl' if no indels are allowed and " +
                "'vj' otherwise.")
cli.c(longOpt: "compress", "Compress output sample files.")


def opt = cli.parse(args)

if (opt == null) {
    System.exit(2)
}

if (opt.h) {
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

def outputFilePrefix = opt.arguments()[-1],
    backgroundSample = opt.b,
    compress = (boolean) opt.c,
    optSearchScope = (opt.o ?: DEFAULT_SEARCH_SCOPE).split(",")

if (optSearchScope.length != 3 || optSearchScope.any { !it.isInteger() || it.toInteger() < 0 }) {
    println "[ERROR] Bad search scope $optSearchScope"
    System.exit(2)
}
def searchScope = optSearchScope.collect { it.toInteger() } as int[]

def grouping = opt.g ?: (searchScope[1] == 0 ? "vjl" : "vj")
grouping = grouping.toLowerCase()
def groupingFactory = grouping == "vjl" ? new VJLClonotypeGroupingFactory() :
        (grouping == "vj" ? new VJClonotypeGroupingFactory() : new DummyClonotypeGroupingFactory())

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Prepare samples

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) : // todo: if no control, perhaps store samples?
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples prepared"

// Compute control degree statistics

def bgDegreeStatCalc = new DegreeStatisticsCalculator(searchScope[0],
        searchScope[1], searchScope[2], groupingFactory)

if (backgroundSample) {
    // Load control sample
    println "[${new Date()} $scriptName] Loading control sample"
    bgDegreeStatCalc.inititalize(SampleFileConnection.load((String) backgroundSample))
} else {
    // Create control by aggregating all samples on the fly
    println "[${new Date()} $scriptName] No control sample provided. Creating pooled sample to be used as control"
    def sampleAggr = new SampleAggregator(sampleCollection,
            new StoringClonotypeAggregatorFactory(),
            OverlapType.Strict)
    bgDegreeStatCalc.inititalize(new PooledSample(sampleAggr))
}

// Compute degree statistics and write output

def sw = new SampleWriter(compress)

sampleCollection.eachWithIndex { sample, ind ->
    def sampleId = sample.sampleMetadata.sampleId
    println "[${new Date()} $scriptName] Computing degree statistics for $sampleId.."

    // Initialize sample degree statistics
    def sampleDegreeStatCalc = new DegreeStatisticsCalculator(searchScope[0],
            searchScope[1], searchScope[2], groupingFactory)
    sampleDegreeStatCalc.inititalize(sample)

    // Create annotator and annotate sample with both background and observed degree statistics
    def annotator = new DegreeStatisticsAnnotator(sampleDegreeStatCalc, bgDegreeStatCalc)
    annotator.annotate(sample)

    // print output
    sw.writeConventional(sample, outputFilePrefix)
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "degstat")

println "[${new Date()} $scriptName] Finished"