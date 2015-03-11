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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.intersection.mds.DiscreteFactorClusterStats
import com.antigenomics.vdjtools.sample.metadata.MetadataTable
import com.antigenomics.vdjtools.util.RUtil

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath

def MEASURE_DEFAULT = "F2", I_TYPE_DEFAULT = "aa"
def cli = new CliBuilder(usage: "BatchIntersectPairPlot [options] input_prefix [output_prefix]\n" +
        "input_prefix should be the same as the output_prefix" +
        "specified during BatchIntersectPair execution")
cli.h("display help message")
cli.m(longOpt: "measure", argName: "string", args: 1,
        "Distance measure to use, allowed values are ${IntersectMetric.allowedNames}. " +
                "[default = $MEASURE_DEFAULT]")
cli.p(longOpt: "plot", "[plotting] Unused.")
cli.f(longOpt: "factor", argName: "string", args: 1,
        "Column name, as in metadata. Factor used to color the plot. [default = no factor]")
cli.n(longOpt: "num-factor", "Factor will be treated as numeric value and gradient plot coloring will be used. " +
        "Factor can still contain non-numeric values which will be converted to NA (grey). " +
        "Should contain at least one numeric value. [default = off]")
cli.l(longOpt: "label", argName: "string", args: 1,
        "Column name, as in metadata. Row values will be used as sample labels. [default = sample_id]")
cli.k(longOpt: "hcl-cutoff", argName: "int, >0", args: 1,
        "Number of clusters to cut the dendrogram into. No output is generated for values <1. [default = 3]")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule, as used in BatchIntersectPairPlot." +
                "Allowed values: $IntersectionType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")

def opt = cli.parse(args)

if (opt == null)
    System.exit(0)

if (opt.h || opt.arguments().size() < 1) {
    cli.usage()
    System.exit(0)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = IntersectionType.getByShortName(iName)

if (!intersectionType) {
    println "[ERROR] Bad intersection type specified ($iName). " +
            "Allowed values are: $IntersectionType.allowedNames"
    System.exit(-1)
}

def inputPrefix = opt.arguments()[0],
    inputFileName = formOutputPath(inputPrefix, "intersect", "batch", intersectionType.shortName)

if (!new File(inputFileName).exists()) {
    println "[ERROR] Input file $inputFileName not found"
    System.exit(-1)
}

intersectionType = intersectionType.shortName

def outputPrefix = opt.arguments().size() > 1 ? opt.arguments()[1] : inputPrefix,
    sampleId = MetadataTable.SAMPLE_ID_COLUMN.toUpperCase(), factorName = opt.f, numFactor = opt.n,
    measureName = (opt.m ?: MEASURE_DEFAULT).toUpperCase(), labelName = (opt.l ?: MetadataTable.SAMPLE_ID_COLUMN).toUpperCase(),
    hcFileName = formOutputPath(outputPrefix, "hc", intersectionType, measureName, "pdf"),
    mdsFileName = formOutputPath(outputPrefix, "mds", intersectionType, measureName, "pdf"),
    k = (opt.k ?: 3),
    clustFileName = formOutputPath(outputPrefix, "hc", "cut$k", intersectionType, measureName),
    coordFileName = formOutputPath(outputPrefix, "mds", "coords", intersectionType, measureName)

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
// we do it manually here as we have no access to metadata :(
boolean specifiedFactor = factorCol1Ind.toInteger() > 0
if (numFactor && specifiedFactor) {
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

//
// Cluster & plot
//

println "[${new Date()} $scriptName] Clustering and plotting data"

RUtil.execute("batch_intersect_pair_clust.r",
        inputFileName,
        idCol1Ind, idCol2Ind,
        measureColInd, IntersectMetric.getByShortName(measureName).normalization.id.toString(),
        factorCol1Ind, factorCol2Ind,
        labelCol1Ind, labelCol2Ind,
        factorNameOrig ?: "NA", numFactor ? "TRUE" : "FALSE",
        hcFileName, mdsFileName,
        k.toString(), clustFileName, coordFileName
)

//
// Permutation testing
//

if (specifiedFactor) {
    if (numFactor) {
        // todo: finish with distance correlation
    } else {
        println "[${new Date()} $scriptName] Running permutation testing for factor ~ cluster dependence"
        def permsOutputPath = formOutputPath(outputPrefix, "perms", intersectionType, measureName)
        def summary = new DiscreteFactorClusterStats(coordFileName).performPermutations(10000)
        if (summary) {
            DiscreteFactorClusterStats.writeSummary(summary, permsOutputPath)
            RUtil.execute("batch_intersect_pair_perm_f.r",
                    permsOutputPath,
                    toPlotPath(permsOutputPath)
            )
            new File(permsOutputPath).delete()
        } else {
            println "[${new Date()} $scriptName] No way - less than 2 factor levels are present."
        }
    }
}

println "[${new Date()} $scriptName] Finished"
