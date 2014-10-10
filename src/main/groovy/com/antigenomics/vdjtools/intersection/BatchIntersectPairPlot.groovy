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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.util.RUtil

def cli = new CliBuilder(usage: "BatchIntersectPairPlot [options] input_file output_prefix")
cli.h("display help message")
cli.m(longOpt: "measure", argName: "string", args: 1,
        "Distance measure to use, vJSD (variable segment Jensen-Shannon distance), " +
                "F (frequency), D (diversity) or R (correlation). [default = F]")
cli.f(longOpt: "factor", argName: "string", args: 1,
        "Column name, as in metadata. Factor used to color the plot. [default = no factor]")
cli.n(longOpt: "num-factor", "Factor will be treated as numeric value and gradient plot coloring will be used. " +
        "Factor can still contain non-numeric values which will be converted to NA (grey). " +
        "Should contain at least one numeric value. [default = off]")
cli.l(longOpt: "label", argName: "string", args: 1,
        "Column name, as in metadata. Row values will be used as sample labels. [default = sample_id]")
cli.k(longOpt: "hcl-cutoff", argName: "int, >0", args: 1,
        "Number of clusters to cut the dendrogram into. No output is generated for values <1. [default = 3]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

def inputFileName = opt.arguments()[0],
    sampleId = "sample_id".toUpperCase(), factorName = opt.f, numFactor = opt.n,
    measureName = (opt.m ?: "F").toUpperCase(), labelName = (opt.l ?: "sample_id").toUpperCase(),
    hcFileName = opt.arguments()[1] + ".batch_intersect_hc.pdf",
    mdsFileName = opt.arguments()[1] + ".batch_intersect_mds.pdf",
    k = (opt.k ?: 3),
    clustFileName = opt.arguments()[1] + ".batch_intersect_clusters${k}.txt"

def factorNameOrig = null
if (factorName) {
    factorNameOrig = factorName
    factorName = factorName.toUpperCase()
}

// Read header

println "[${new Date()} $scriptName] Reading data header"

def header = []
new File(inputFileName).withReader { reader ->
    header = reader.readLine().split("\t").collect { it.toUpperCase() }
}

// Match column indices

def idCol1Ind = (header.findIndexOf { it.contains("1_$sampleId") } + 1).toString(),
    idCol2Ind = (header.findIndexOf { it.contains("2_$sampleId") } + 1).toString(),
    measureColInd = (header.findIndexOf { it.equals(measureName) } + 1).toString(),
    factorCol1Ind = ((factorName ? header.findIndexOf { it.contains("1_$factorName") } : -1) + 1).toString(),
    factorCol2Ind = ((factorName ? header.findIndexOf { it.contains("2_$factorName") } : -1) + 1).toString(),
    labelCol1Ind = (header.findIndexOf { it.contains("1_$labelName") } + 1).toString(),
    labelCol2Ind = (header.findIndexOf { it.contains("2_$labelName") } + 1).toString()

if (measureColInd.toInteger() < 1) {
    println "[ERROR] Measure column ($measureName) is absent. Terminating"
    System.exit(-1)
}

// Check if we can map factor to gradient scale
if (numFactor && factorCol1Ind.toInteger() > 0) {
    int fcol1 = factorCol1Ind.toInteger() - 1, fcol2 = factorCol2Ind.toInteger() - 1
    def fValues = new HashSet<Double>()

    new File(inputFileName).withReader { reader ->
        header = reader.readLine().split("\t")
        def line
        while ((line = reader.readLine()) != null) {
            def splitLine = line.split("\t")
            if (splitLine[fcol1].isDouble())
                fValues.add(splitLine[fcol1].toDouble())
            if (splitLine[fcol2].isDouble())
                fValues.add(splitLine[fcol2].toDouble())
        }
    }

    if (fValues.size() < 3) {
        println "[WARNING] Numeric factor specified, while number of unique numeric factor values < 3. Switching it off.."
        numFactor = false
    }
}
// Plot

println "[${new Date()} $scriptName] Plotting data"

RUtil.execute("batch_intersect_pair_clust.r",
        inputFileName,
        idCol1Ind, idCol2Ind,
        measureColInd,
        factorCol1Ind, factorCol2Ind,
        labelCol1Ind, labelCol2Ind,
        factorNameOrig ?: "NA", numFactor ? "TRUE" : "FALSE",
        hcFileName, mdsFileName,
        clustFileName, k.toString()
)


