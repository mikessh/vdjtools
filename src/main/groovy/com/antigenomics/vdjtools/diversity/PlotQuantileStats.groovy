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
 */

package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.RUtil

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath

def TOP_DEFAULT = "5", TOP_MAX = 10
def cli = new CliBuilder(usage: "PlotQuantileStats [options] input_name output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.t(longOpt: "top", args: 1, "Number of top clonotypes to present on the histogram. " +
        "Values > $TOP_MAX are not allowed, as they would make the plot unreadable. [default = $TOP_DEFAULT]")

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

//
// Calc stats
//

println "[${new Date()} $scriptName] Computing quantile stats"

def quantileStats = new QuantileStats(sample)

//
// Write output table
//

println "[${new Date()} $scriptName] Writing output"

def outputFileName = formOutputPath(outputFilePrefix, "qstat")

new File(outputFileName).withPrintWriter { pw ->
    pw.println(QuantileStats.HEADER)
    pw.println(quantileStats)
    (0..<top).each {
        def clonotype = sample[it]
        pw.println("top\t${clonotype.cdr3aa}\t${clonotype.freq}")
    }
}

//
// Create sunburst plot
//

println "[${new Date()} $scriptName] Plotting"

RUtil.execute("quantile_stats.r",
        outputFileName,
        toPlotPath(outputFileName)
)

println "[${new Date()} $scriptName] Finished"



