package com.antigenomics.vdjtools.diversity

import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.intersection.IntersectionType
import com.antigenomics.vdjtools.sample.SampleCollection
import com.antigenomics.vdjtools.util.ExecUtil
import com.antigenomics.vdjtools.util.RUtil

/**
 * Created by mikesh on 10/30/14.
 */

def I_TYPE_DEFAULT = "strict"
def cli = new CliBuilder(usage: "BuildFrequencyTable [options] input_name output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Comma-separated list of intersection types to apply. " +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.p(longOpt: "plot", "Plots rarefaction curves. " +
        "(R installation with ggplot2 and reshape packages is reuqired).")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S), plot = opt.p,
    outputPrefix = opt.arguments()[1]

ExecUtil.ensureDir(outputPrefix)

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

new File(outputPrefix + ".freqtable.txt").withPrintWriter { pw ->
    pw.println(FrequencyTable.BinInfo.HEADER)
    frequencyTable.bins.each {
        pw.println(it)
    }
}

if (plot) {
    println "[${new Date()} $scriptName] Plotting data"

    RUtil.execute("freqtable_plot.r",
            outputPrefix + ".freqtable.txt",
            outputPrefix + ".freqtable.pdf"
    )
}

println "[${new Date()} $scriptName] Finished"