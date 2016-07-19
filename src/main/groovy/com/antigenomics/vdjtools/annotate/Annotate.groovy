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

import com.antigenomics.vdjtools.io.SampleWriter
import com.antigenomics.vdjtools.sample.SampleCollection

def DEFAULT_BASE_ANNOTS = ["cdr3length", "ndnsize", "insertsize"].join(","),
    DEFAULT_AAPROP_ANNOTS = ["hydropathy", "charge", "polarity", "strength", "high_contact"].join(",")

def cli = new CliBuilder(usage: "Annotate [options] " +
        "[sample1 sample2 ... if not -m] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata.")
cli.b(longOpt: "base", argName: "param1,param2,...", args: 1,
        "Comma-separated list of basic clonotype features to calculate and append to " +
                "resulting clonotype tables. " +
                "Allowed values: ${BaseClonotypeAnnotator.ALLOWED_NAMES} (case-sensitive!). " +
                "[default=$DEFAULT_BASE_ANNOTS]")
cli.a(longOpt: "aaprop", argName: "param1,param2,...", args: 1,
        "Comma-separated list of amino acid properties. Amino acid property value sum will be " +
                "calculated for CDR3 sequence of each coding clonotype and appended to " +
                "resulting clonotype tables (for non-coding clonotypes corresponding fields will be left blank). " +
                "Allowed values: ${AminoAcidPropertyClonotypeAnnotator.ALLOWED_NAMES}. " +
                "[default=$DEFAULT_AAPROP_ANNOTS]")
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
    baseAnnot = (opt.b ?: DEFAULT_BASE_ANNOTS).split(",").collect { it.toLowerCase() },
    aapropAnnot = (opt.a ?: DEFAULT_AAPROP_ANNOTS).split(",").collect { it.toLowerCase() },
    compress = (boolean) opt.c

def scriptName = getClass().canonicalName.split("\\.")[-1]

def badBaseField = baseAnnot.findAll { !BaseClonotypeAnnotator.ALLOWED_NAMES.contains(it) }
badBaseField.addAll(aapropAnnot.findAll { !AminoAcidPropertyClonotypeAnnotator.ALLOWED_NAMES.contains(it) })

if (!badBaseField.empty) {
    println "The following annotation properties are specified incorrectly: $badBaseField"
    System.exit(2)
}

//
// Batch load samples
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples prepared"

//
// Create annotators
//

def baseAnnotators = baseAnnot.collect { new BaseClonotypeAnnotator(it) },
    aapropAnnotators = aapropAnnot.collect { new AminoAcidPropertyClonotypeAnnotator(it) }

def sampleAnnotator = new SampleAnnotator([baseAnnotators, aapropAnnotators].flatten())

//
// Iterate over samples and annotate
//
def sw = new SampleWriter(compress)

sampleCollection.eachWithIndex { sample, ind ->
    def sampleId = sample.sampleMetadata.sampleId
    println "[${new Date()} $scriptName] Annotating $sampleId.."

    sampleAnnotator.annotate(sample)

    // print output
    sw.writeConventional(sample, outputFilePrefix)
}

sampleCollection.metadataTable.storeWithOutput(outputFilePrefix, compress,
        "annot:${baseAnnot.join(",")}:${aapropAnnot.join(",")}")

println "[${new Date()} $scriptName] Finished"
