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

package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Edge
import com.antigenomics.vdjtools.segment.SegmentUtil
import com.antigenomics.vdjtools.segment.VSegmentTable

def FREQ_THRESHOLD = "0.4", SPEC_THRESHOLD = "3"
def cli = new CliBuilder(usage: "IgBlastStat [options] igblast_output_level2 output_prefix")
cli.h("display help message")
cli.S(longOpt: "species", argName: "string",
        "Species for which partitioning info on IG regions (FWs and CDRs) will be loaded. " +
                "Possible values: human [default], mouse, rabbit and rat.")
cli._(longOpt: "allele-freq", argName: "[0, 1]", "Frequency threshold, used together with spectratype threshold. " +
        "Mutations with higher frequency are considered as allele candidates. [default=$FREQ_THRESHOLD]")
cli._(longOpt: "allele-spectra", argName: "int", "Spectratype threshold, used together with frequency threshold. " +
        "At least \$allele-spectra clonotypes with distinct CDR3 lengths must contain this mutation " +
        "for it to be considered as allele. [default=$SPEC_THRESHOLD]")
def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

def species = (opt.S ?: "human").toLowerCase(),
    freqThreshold = (opt.'allele-freq' ?: FREQ_THRESHOLD).toDouble(),
    spectraThreshold = (opt.'allele-spectra' ?: SPEC_THRESHOLD).toInteger()

String inputFileNameL2 = opt.arguments()[0],
       outputPrefix = opt.arguments()[1]

if (!new File(inputFileNameL2).exists()) {
    println "[ERROR] Input file does not exist: $inputFileNameL2"
}

def vSegmentTable = new VSegmentTable(species)

//
// Read and parse clonotypes
//

println "[${new Date()} $scriptName] Loading L2 clonotypes"

def clonotypes = new LinkedList<Clonotype>()
new File(inputFileNameL2).eachLine { line ->
    if (!line.startsWith("#")) {
        clonotypes.add(Clonotype.parseIgBlastClonotype(line))
    }
}

def clonotypeMap = new ClonotypeMap(vSegmentTable, clonotypes)

println "[${new Date()} $scriptName] Deducing alleles"

clonotypeMap.deduceAlleles(freqThreshold, spectraThreshold)

//
// Iterate through clonotypes with same CDR3nt and build hypermutation links
//

println "[${new Date()} $scriptName] Calculating SHM edges for germline"

def shmGraph = new ShmGraphBuilder(clonotypeMap).buildGraph()

//
// CDR3 hypermutations
//

println "[${new Date()} $scriptName] Calculating SHM edges for CDR3"

def cdr3Graph = new Cdr3GraphBuilder(clonotypeMap).buildGraph()


def of = new File(outputPrefix).absoluteFile
if (of.parentFile != null)
    of.parentFile.mkdirs()

//
// Generate tables
//

def alleles = clonotypeMap.clonotypes.collect { it.alleles }.flatten(),
    shms = clonotypeMap.clonotypes.collect { it.shms }.flatten(),
    shmsEmerged = [
            shmGraph.filteredShms,
            cdr3Graph.filteredShms
    ].flatten()

// RS table

println "[${new Date()} $scriptName] Calculating Silent:Replacement tables"

def allelesRsTable = new RSTable(true, vSegmentTable)
allelesRsTable.addAll(alleles)

def shmRsTable = new RSTable(true, vSegmentTable)
shmRsTable.addAll(shms)

def emergedShmRsTable = new RSTable(true, vSegmentTable)
emergedShmRsTable.addAll(shmsEmerged)

// Mutation motifs

println "[${new Date()} $scriptName] Calculating SHM hotspot motifs"

def allelesPwm = new MotifPwm(), shmPwm = new MotifPwm(), emergedShmPwm = new MotifPwm()
allelesPwm.addAll(alleles)
shmPwm.addAll(shms)
emergedShmPwm.addAll(shmsEmerged)

//
// Write output
//

println "[${new Date()} $scriptName] Writing output"

new File(outputPrefix + ".mutations.rs.txt").withPrintWriter { pw ->
    pw.println("#silent:replacement\t" + SegmentUtil.HEADER)
    pw.println("alleles\t" + allelesRsTable.summaryRs.collect().join("\t"))
    pw.println("shms\t" + shmRsTable.summaryRs.collect().join("\t"))
    pw.println("shms_emerged\t" + emergedShmRsTable.summaryRs.collect().join("\t"))
}

new File(outputPrefix + ".mutations.cov.txt").withPrintWriter { pw ->
    pw.println("#coverage\t" + SegmentUtil.HEADER)
    pw.println("alleles\t" + allelesRsTable.summaryCoverage.collect().join("\t"))
    pw.println("shms\t" + shmRsTable.summaryCoverage.collect().join("\t"))
    pw.println("shms_emerged\t" + emergedShmRsTable.summaryCoverage.collect().join("\t"))
}

new File(outputPrefix + ".mutations.pwm.txt").withPrintWriter { pw ->
    pw.println("#alleles")
    pw.println(allelesPwm)
    pw.println("\n#shms")
    pw.println(shmPwm)
    pw.println("\n#emerged_shms")
    pw.println(emergedShmPwm)
}

// Cytoscape files

new File(outputPrefix + ".nodes.txt").withPrintWriter { pw ->
    pw.println(Clonotype.NODE_HEADER)
    clonotypeMap.clonotypes.each {
        pw.println(it.nodeString())
    }
}

new File(outputPrefix + ".edges.txt").withPrintWriter { pw ->
    pw.println(Edge.EDGE_HEADER)

    shmGraph.edges.each {
        pw.println(it.edgeString())
    }

    cdr3Graph.edges.each {
        pw.println(it.edgeString())
    }
}

new File(outputPrefix + ".net.txt").withPrintWriter { pw ->
    pw.println(Edge.NET_HEADER)

    shmGraph.edges.each {
        pw.println(it.netString())
    }

    cdr3Graph.edges.each {
        pw.println(it.netString())
    }
}

println "[${new Date()} $scriptName] Finished"