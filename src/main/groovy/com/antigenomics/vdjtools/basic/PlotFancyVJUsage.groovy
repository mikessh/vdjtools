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

def cli = new CliBuilder(usage: "PlotFancyVJUsage [options] input_name output_prefix")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")
cli.u(longOpt: "unweighted", "Will count each clonotype only once, apart from conventional frequency-weighted histogram.")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S),
    unweighted = opt.u,
    outputPrefix = opt.arguments()[1]

ExecUtil.ensureDir(outputPrefix)

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Read the sample
//

println "[${new Date()} $scriptName] Reading sample"

def sampleCollection = new SampleCollection([opt.arguments()[0]], software)

def sampleId = sampleCollection.metadataTable.getRow(0).sampleId

// Calculate segment usage
def segmentUsage = new SegmentUsage(sampleCollection, unweighted)

// Output and plotting
println "[${new Date()} $scriptName] Writing output"

new File(outputPrefix + ".fancyvj.txt").withPrintWriter { pw ->
    pw.println(".\t" + segmentUsage.vUsageHeader().collect().join("\t"))
    def vjMatrix = segmentUsage.vjUsageMatrix(sampleId)
    vjMatrix.eachWithIndex { double[] vVectorByJ, int i ->
        pw.println(segmentUsage.jUsageHeader()[i] + "\t" + vVectorByJ.collect().join("\t"))
    }
}

println "[${new Date()} $scriptName] Plotting data (be patient, complex graphics)"

RUtil.execute("table_circ.r",
        outputPrefix + ".fancyvj.txt", outputPrefix + ".fancyvj.pdf"
)


println "[${new Date()} $scriptName] Finished"