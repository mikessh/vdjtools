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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil
import com.antigenomics.vdjtools.util.RUtil


def cli = new CliBuilder(usage: "PlotFancySpectratype [options] input_name output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.t(longOpt: "top", args: 1, "Number of top clonotypes to present on the histogram. " +
        "Values > 20 are not allowed, as they would make the plot legend unreadable. default = 20")


def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S),
    outputPrefix = opt.arguments()[1],
    top = opt.t ?: 20

if (top > 20) {
    println "[ERROR] Specified number of top clonotypes should not exceed 20"
    System.exit(-1)
}

ExecUtil.ensureDir(outputPrefix)

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Read the sample
//

println "[${new Date()} $scriptName] Reading sample"

def sampleCollection = new SampleCollection([opt.arguments()[0]], software, false)

def sample = sampleCollection[0]
sample.renormalize() // ensure normalization as we're going to subtract top clonotype frequencies

// Calculate spectratype
def spectratype = new Spectratype(aminoAcid, false)

spectratype.addAll(sample)

def spectratypeHist = spectratype.histogram

// Calculate top clonotypes and subtract their frequencies

def topClonotypes = sample.top(top)
topClonotypes.each {
    def bin = spectratype.bin(it)
    spectratypeHist[bin] -= it.freq
}

// Prepair output table

def spectraMatrix = new double[spectratypeHist.length][top + 1]

for (int i = 0; i < spectratypeHist.length; i++) {
    spectraMatrix[i][0] = spectratypeHist[i]
}

topClonotypes.sort { it.freq }.eachWithIndex { it, ind ->
    def bin = spectratype.bin(it)
    spectraMatrix[bin][ind + 1] = it.freq
}

def table = "Len\tOther\t" + topClonotypes.collect { it.cdr3aa }.join("\t")
for (int i = 0; i < spectratypeHist.length; i++) {
    table += "\n" + spectratype.lengths[i] + "\t" + spectraMatrix[i].collect().join("\t")
}

// Output

println "[${new Date()} $scriptName] Writing output and plotting data"

new File(outputPrefix + ".fancyspectra.txt").withPrintWriter { pw ->
    pw.println(table)
}

RUtil.execute("fancy_spectratype.r",
        table, outputPrefix + ".fancyspectra.pdf"
)

println "[${new Date()} $scriptName] Finished"


