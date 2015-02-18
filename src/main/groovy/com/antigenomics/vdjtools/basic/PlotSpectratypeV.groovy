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

package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.RUtil

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath

def TOP_DEFAULT = "12", TOP_MAX = 12
def cli = new CliBuilder(usage: "PlotSpectratypeV [options] input_name output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.t(longOpt: "top", args: 1, "Number of top V segments to present on the histogram. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot unreadable. [default = $TOP_DEFAULT]")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S),
    outputFilePrefix = opt.arguments()[1],
    top = (opt.t ?: TOP_DEFAULT).toInteger(),
    unweighted = (boolean) opt.u

if (top > TOP_MAX) {
    println "[ERROR] Specified number of top V segments should not exceed 20"
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

def spectratypeV = new SpectratypeV(false, unweighted)

spectratypeV.addAll(sample)

def collapsedSpectratypes = spectratypeV.collapse(top)

// Prepare output table

def spectraMatrix = new double[spectratypeV.span][top + 1]

def otherHistogram = collapsedSpectratypes["other"].getHistogram(false)
for (int i = 0; i < spectratypeV.span; i++) {
    spectraMatrix[i][0] = otherHistogram[i]
}

collapsedSpectratypes.findAll { it.key != "other" }.eachWithIndex { it, ind ->
    def histogram = it.value.getHistogram(false)
    for (int i = 0; i < spectratypeV.span; i++) {
        spectraMatrix[i][top - ind] = histogram[i]
    }
}

def table = "Len\tOther\t" + collapsedSpectratypes.findAll { it.key != "other" }.collect { it.key }.reverse().join("\t")
for (int i = 0; i < spectratypeV.span; i++) {
    table += "\n" + spectratypeV.lengths[i] + "\t" + spectraMatrix[i].collect().join("\t")
}

// Output

println "[${new Date()} $scriptName] Writing output and plotting data"

def outputFileName = formOutputPath(outputFilePrefix, "spectraV", (unweighted ? "unwt" : "wt"))

new File(outputFileName).withPrintWriter { pw ->
    pw.println(table)
}

RUtil.execute("fancy_spectratype.r",
        outputFileName, toPlotPath(outputFileName), "Variable segment", RUtil.logical(false)
)

println "[${new Date()} $scriptName] Finished"