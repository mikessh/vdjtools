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

package com.antigenomics.vdjtools.overlap

import com.antigenomics.vdjtools.overlap.permutations.DiscreteFactorClusterStats

import static com.antigenomics.vdjtools.util.ExecUtil.formOutputPath
import static com.antigenomics.vdjtools.util.ExecUtil.toPlotPath
import static com.antigenomics.vdjtools.util.RUtil.execute

def MEASURE_DEFAULT = "F", I_TYPE_DEFAULT = "aa"

def cli = new CliBuilder(usage: "TestClusters [options] input_prefix [output_prefix]\n" +
        "NOTE: input_prefix should be equal to output_prefix specified for" +
        "ClusterSamples execution, -i, -e parameters should also match.")

cli.h("display help message")
cli.e(longOpt: "measure", argName: "string", args: 1,
        "Distance measure to use, allowed values are ${OverlapMetric.allowedNames}. " +
                "[default = $MEASURE_DEFAULT]")
cli.i(longOpt: "intersect-type", argName: "string", args: 1,
        "Intersection rule, as used in CalcPairwiseDistances." +
                "Allowed values: $OverlapType.allowedNames. " +
                "Will use '$I_TYPE_DEFAULT' by default.")
cli.n(longOpt: "num-factor", "Treat factor as numeric")
cli._(longOpt: "plot-type", argName: "<pdf|png>", args: 1, "Plot output format [default=pdf]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(0)

if (opt.h || opt.arguments().size() < 1) {
    cli.usage()
    System.exit(0)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

def iName = opt.i ?: I_TYPE_DEFAULT
def intersectionType = OverlapType.getByShortName(iName)

if (!intersectionType) {
    println "[ERROR] Bad overlap type specified ($iName). " +
            "Allowed values are: $OverlapType.allowedNames"
    System.exit(-1)
}

intersectionType = intersectionType.shortName

def measureName = (opt.e ?: MEASURE_DEFAULT).toUpperCase(),
    numFactor = (boolean) opt.n,
    inputPrefix = opt.arguments()[0],
    mdsFileName = formOutputPath(inputPrefix, "mds", intersectionType, measureName),
    outputPrefix = opt.arguments().size() > 1 ? opt.arguments()[1] : inputPrefix,
    plotType = (opt.'plot-type' ?: "pdf").toString()

//
// Permutation testing
//

if (numFactor) {
    // todo: finish with distance correlation
} else {
    println "[${new Date()} $scriptName] Running permutation testing for factor ~ cluster dependence"
    def permsOutputPath = formOutputPath(outputPrefix, "perms", intersectionType, measureName)
    def summary = new DiscreteFactorClusterStats(mdsFileName).performPermutations(10000)
    if (summary) {
        DiscreteFactorClusterStats.writeSummary(summary, permsOutputPath)
        execute("cluster_permutations_plot.r",
                permsOutputPath,
                toPlotPath(permsOutputPath, plotType)
        )
        new File(permsOutputPath).delete()
    } else {
        println "[${new Date()} $scriptName] No way - less than 2 factor levels are present."
    }
}

println "[${new Date()} $scriptName] Finished"
