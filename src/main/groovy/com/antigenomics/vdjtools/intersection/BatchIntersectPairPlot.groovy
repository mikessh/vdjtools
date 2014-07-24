package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.RUtil

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

def cli = new CliBuilder(usage: "BatchIntersectPairPlot [options] input_file output_prefix")
cli.h("display help message")
cli.m(longOpt: "measure", argName: "string", args: 1,
        "Distance measure to use, F (frequency), D (diversity) or R (correlation). [default = F]")
cli.f(longOpt: "factor", argName: "string", args: 1,
        "Column name, as in metadata. Factor used to color the plot. " +
                "Can contain non-numeric values which will be converted to NA. " +
                "Should contain at least one numeric value. [default = no factor]")
cli.l(longOpt: "label", argName: "string", args: 1,
        "Column name, as in metadata. Row values will be used as sample labels. [default = sample_id]")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

/*file_in           = "intersection2_aa.txt"
id_col1_index     = 1#"X.1_sample_id"
id_col2_index     = 9#"X2_sample_id"
measure_col_index = 22#"freq12"
factor_col1_index = 4#"X1_age"
factor_col2_index = 12#"X2_age"
lbl_col1_index    = 4  #
lbl_col2_index    = 12 #*/

def inputFileName = opt.arguments()[0],
    sampleId = "sample_id", factorName = opt.f,
    measureName = (opt.m ?: "F"), labelName = (opt.l ?: "sample_id"),
    hcFileName = opt.arguments()[1] + "_hc.pdf", mdsFileName = opt.arguments()[1] +  "_mds.pdf"

// Read header

println "[${new Date()} $scriptName] Reading data header"

def header = []
new File(inputFileName).withReader { reader ->
    header = reader.readLine().split("\t")
}

// Match column indices
def idCol1Ind = (header.findIndexOf { it.contains("1_$sampleId") } + 1).toString(),
    idCol2Ind = (header.findIndexOf { it.contains("2_$sampleId") } + 1).toString(),
    measureColInd = (header.findIndexOf { it.contains(measureName) } + 1).toString(),
    factorCol1Ind = ((factorName ? header.findIndexOf { it.contains("1_$factorName") } : -1) + 1).toString(),
    factorCol2Ind = ((factorName ? header.findIndexOf { it.contains("2_$factorName") } : -1) + 1).toString(),
    labelCol1Ind = (header.findIndexOf { it.contains("1_$labelName") } + 1).toString(),
    labelCol2Ind = (header.findIndexOf { it.contains("2_$labelName") } + 1).toString()

// Plot
println "[${new Date()} $scriptName] Plotting data"

RUtil.execute("batch_intersect_pair_clust.r",
        inputFileName,
        idCol1Ind, idCol2Ind,
        measureColInd,
        factorCol1Ind, factorCol2Ind,
        labelCol1Ind, labelCol2Ind,
        hcFileName, mdsFileName
)


