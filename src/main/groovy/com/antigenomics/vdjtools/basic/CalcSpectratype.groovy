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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.sample.metadata.MetadataTable

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath

def cli = new CliBuilder(usage: "CalcSpectratype [options] " +
        "[sample1 sample2 sample3 ... if -m is not specified] output_prefix")
cli.h("display help message")
cli.m(longOpt: "metadata", argName: "filename", args: 1,
        "Metadata file. First and second columns should contain file name and sample id. " +
                "Header is mandatory and will be used to assign column names for metadata." +
                "If column named 'time' is present, it will be used to specify time point sequence.")
cli.a(longOpt: "amino-acid", "Will use amino-acid CDR3 sequence lengths instead of nucleotide.")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")

def opt = cli.parse(args)

if (opt == null || opt.arguments().size() == 0)
    System.exit(-1)

if (opt.h) {
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

def outputFilePrefix = opt.arguments()[-1],
    aminoAcid = (boolean) opt.a, unweighted = (boolean) opt.u


def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Batch load all samples (lazy)
//

println "[${new Date()} $scriptName] Reading samples"

def sampleCollection = metadataFileName ?
        new SampleCollection((String) metadataFileName) :
        new SampleCollection(opt.arguments()[0..-2])

println "[${new Date()} $scriptName] ${sampleCollection.size()} samples loaded"

//
// Compute and output diversity measures, spectratype, etc
//

new File(formOutputPath(outputFilePrefix, "spectratype", (aminoAcid ? "aa" : "nt"), (unweighted ? "unwt" : "wt"))).withPrintWriter { pwSpectra ->
    new File(formOutputPath(outputFilePrefix, "spectratype.insert", (unweighted ? "unwt" : "wt"))).withPrintWriter { pwIns ->
        new File(formOutputPath(outputFilePrefix, "spectratype.ndn", (unweighted ? "unwt" : "wt"))).withPrintWriter { pwNdn ->
            def spectratype = new Spectratype(aminoAcid, unweighted),
                insertHist = new InsertHist(unweighted),
                ndnHist = new NdnHist(unweighted)

            def header = "#$MetadataTable.SAMPLE_ID_COLUMN\t" + sampleCollection.metadataTable.columnHeader + "\t"

            pwSpectra.println(header + spectratype.HEADER)
            pwIns.println(header + insertHist.HEADER)
            pwNdn.println(header + ndnHist.HEADER)

            def sampleCounter = 0

            sampleCollection.each { Sample sample ->
                spectratype.addAll(sample)
                insertHist.addAll(sample)
                ndnHist.addAll(sample)

                println "[${new Date()} $scriptName] ${++sampleCounter} samples processed"

                pwSpectra.println([sample.sampleMetadata.sampleId, sample.sampleMetadata, spectratype].join("\t"))
                pwIns.println([sample.sampleMetadata.sampleId, sample.sampleMetadata, insertHist].join("\t"))
                pwNdn.println([sample.sampleMetadata.sampleId, sample.sampleMetadata, ndnHist].join("\t"))

                spectratype.clear()
                insertHist.clear()
                ndnHist.clear()
            }
        }
    }
}

println "[${new Date()} $scriptName] Finished"
