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
 * Last modified on 2.11.2014 by mikesh
 */



package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.SampleCollection
import static com.antigenomics.vdjtools.util.ExecUtil.*
import com.antigenomics.vdjtools.util.RUtil

def I_TYPE_DEFAULT = "strict"
def cli = new CliBuilder(usage: "BuildFrequencyTable [options] input_name output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.p(longOpt: "plot", "Plots log-log plot for frequency distribution. " +
        "(R installation with ggplot2 and reshape packages is required).")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), plot = opt.p,
    outputPrefix = opt.arguments()[1]

def scriptName = getClass().canonicalName.split("\\.")[-1]

// Select intersection type

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = IntersectionType.byName(iName)

if (!intersectionType) {
    println "[ERROR] Bad intersection type specified ($iName). " +
            "Allowed values are: $IntersectionType.allowedNames"
    System.exit(-1)
}

//
// Read the sample
//

println "[${new Date()} $scriptName] Reading sample"

def sampleCollection = new SampleCollection([opt.arguments()[0]], software)

def sample = sampleCollection[0]

println "[${new Date()} $scriptName] Building frequency table"

def frequencyTable = new FrequencyTable(sample, intersectionType)

println "[${new Date()} $scriptName] Writing output"

def outputTablePath = formOutputPath(outputPrefix, "freqtable")

new File(outputTablePath).withPrintWriter { pw ->
    pw.println(FrequencyTable.BinInfo.HEADER)
    frequencyTable.bins.each {
        pw.println(it)
    }
}

if (plot) {
    println "[${new Date()} $scriptName] Plotting data"

    RUtil.execute("freqtable_plot.r",
            outputTablePath,
            toPlotPath(outputTablePath)
    )
}

println "[${new Date()} $scriptName] Finished"