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
 * Last modified on 18.10.2014 by mikesh
 */

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.SampleCollection
import static com.antigenomics.vdjtools.util.ExecUtil.*
import com.antigenomics.vdjtools.util.RUtil

def TOP_DEFAULT = "10", TOP_MAX = 20
def cli = new CliBuilder(usage: "PlotFancySpectratype [options] input_name output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.t(longOpt: "top", args: 1, "Number of top clonotypes to present on the histogram. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot legend unreadable. [default = $TOP_DEFAULT]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S),
    outputFilePrefix = opt.arguments()[1],
    top = (opt.t ?: TOP_DEFAULT).toInteger()

if (top > TOP_MAX) {
    println "[ERROR] Specified number of top clonotypes should not exceed $TOP_MAX"
    System.exit(-1)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Read the sample
//

println "[${new Date()} $scriptName] Reading sample"

def sampleCollection = new SampleCollection([opt.arguments()[0]], software)

def sample = sampleCollection[0]


// Calculate spectratype

def spectratype = new Spectratype(false, false)

def topClonotypes = spectratype.addAllFancy(sample, top)

def spectratypeHist = spectratype.histogram


// Prepare output table

def spectraMatrix = new double[spectratype.len][top + 1]

for (int i = 0; i < spectratype.len; i++) {
    spectraMatrix[i][0] = spectratypeHist[i]
}

topClonotypes.eachWithIndex { it, ind ->
    def bin = spectratype.bin(it)
    spectraMatrix[bin][top - ind] = it.freq
}

def table = "Len\tOther\t" + topClonotypes.reverse().collect { it.cdr3aa }.join("\t")
for (int i = 0; i < spectratype.len; i++) {
    table += "\n" + spectratype.lengths[i] + "\t" + spectraMatrix[i].collect().join("\t")
}


// Output

println "[${new Date()} $scriptName] Writing output and plotting data"

def outputFileName = formOutputPath(outputFilePrefix, "fancyspectra")

new File(outputFileName).withPrintWriter { pw ->
    pw.println(table)
}

RUtil.execute("fancy_spectratype.r",
        outputFileName, toPlotPath(outputFileName), "Clonotype", "TRUE"
)

println "[${new Date()} $scriptName] Finished"


