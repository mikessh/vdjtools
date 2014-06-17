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
import com.antigenomics.vdjtools.Mutation
import com.antigenomics.vdjtools.Util
import com.antigenomics.vdjtools.segment.SegmentUtil
import com.antigenomics.vdjtools.segment.VSegmentTable

def FREQ_THRESHOLD = "0.4", SPEC_THRESHOLD = "3", V_FREQ_THRESHOLD = 0.01
def cli = new CliBuilder(usage: "IgBlastStat [options] igblast_output_level2 output_prefix")
cli.h("display help message")
cli.S(longOpt: "species", argName: "string",
        "Species for which partitioning info on IG regions (FWs and CDRs) will be loaded. " +
                "Possible values: human [default], mouse, rabbit and rat.")
cli.d("Deduce CDR3 SHM direction [experimental]")
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
    spectraThreshold = (opt.'allele-spectra' ?: SPEC_THRESHOLD).toInteger(),
    deduceDirection = opt.d

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

def clonotypeMap = new ClonotypeMap(vSegmentTable, inputFileNameL2)

println "[${new Date()} $scriptName] Deducing alleles"

clonotypeMap.deduceAlleles(freqThreshold, V_FREQ_THRESHOLD, spectraThreshold)

//
// Iterate through clonotypes with same CDR3nt and build hypermutation links
//

def cl2clShmMap = new HashMap<String, Map<String, List<Mutation>>>()
def edgeWeightMap = new HashMap<String, Integer>(),
    inDegreeMap = new HashMap<String, Integer>()

/**
 * This procedure puts key2 -> key1 link to mutation network and updates the degrees
 */
def putShm = { String key1, String key2, Mutation shm ->
    def clShmMap = cl2clShmMap[key2]
    if (clShmMap == null)
        cl2clShmMap.put(key2,
                clShmMap = new HashMap<String, List<Mutation>>())
    def shmList = clShmMap[key1]
    if (shmList == null)
        clShmMap.put(key1, shmList = new LinkedList<Mutation>())
    shmList.add(shm.reassignParent(clonotypeMap.getByKey(key1)))

    // Compute incoming degrees
    inDegreeMap.put(key1, (inDegreeMap[key1] ?: 0) + 1)

    // Weight (for visualization, undirected)
    def key = key1 + "\t" + key2
    edgeWeightMap.put(key, (edgeWeightMap[key] ?: 0) + 1)
    key = key2 + "\t" + key1
    edgeWeightMap.put(key, (edgeWeightMap[key] ?: 0) + 1)
}

println "[${new Date()} $scriptName] Creating SHM links for FW1-FW3"

clonotypeMap.clonotypesByCdr3.each { clonotypeList ->
    int maxLevel = 0

    clonotypeList.each { maxLevel = Math.max(maxLevel, 2 * it.shms.size()) }

    final def connectivityMap = new HashMap<String, List<String>>(),
              newConnectivityMap = new HashMap<String, List<String>>()

    def addClonotypePair = { String key1, String key2 ->
        def conList = newConnectivityMap[key1]
        if (conList == null)
            newConnectivityMap.put(key1, conList = new LinkedList<String>())
        conList.add(key2)
    }

    // Pre-compute lists of intersecting SHMs
    def shmIntersections = new Set[clonotypeList.size()][clonotypeList.size()]
    for (int i = 0; i < clonotypeList.size(); i++) {
        def clonotype1 = clonotypeList[i]
        for (int j = i + 1; j < clonotypeList.size(); j++) {
            def clonotype2 = clonotypeList[j]
            def shms1 = clonotype1.shms, shms2 = clonotype2.shms

            def shmCount1 = shms1.size(),
                shmCount2 = shms2.size()

            def shmIntersection = shmCount1 > shmCount2 ?
                    shms1.intersect(shms2) : shms2.intersect(shms1)

            shmIntersections[i][j] = shmIntersection
        }
    }

    // Iteratively scan for differences by 1,2,.. SHMs
    // Check connectivity graph at previous iteration to see if we're not creating cycles
    for (int level = 1; level < maxLevel; level++) {
        for (int i = 0; i < clonotypeList.size(); i++) {
            def clonotype1 = clonotypeList[i]
            for (int j = i + 1; j < clonotypeList.size(); j++) {
                def clonotype2 = clonotypeList[j]

                def shms1 = clonotype1.shms, shms2 = clonotype2.shms
                def shmCount1 = shms1.size(),
                    shmCount2 = shms2.size()
                def shmIntersection = shmIntersections[i][j]
                def shmDelta1 = shmCount1 - shmIntersection.size(),
                    shmDelta2 = shmCount2 - shmIntersection.size()

                def key1 = clonotype1.key, key2 = clonotype2.key

                if (Math.abs(shmDelta1) + Math.abs(shmDelta2) == level) {
                    boolean connected = false

                    def connectivityCheck = new ConnectivityCheck()

                    if (!connectivityCheck.checkNodes(connectivityMap, key1, key2)) {
                        connected = true
                        shms1.each { shm ->
                            if (!shmIntersection.contains(shm))
                                putShm(key1, key2, shm)
                        }
                        shms2.each { shm ->
                            if (!shmIntersection.contains(shm))
                                putShm(key2, key1, shm)
                        }
                    }

                    // update connectivity map
                    if (connected) {
                        addClonotypePair(key1, key2)
                        addClonotypePair(key2, key1)
                    }
                }
            }
        }

        // Merge connectivity map
        newConnectivityMap.each {
            def conList = connectivityMap[it.key]
            if (conList == null)
                connectivityMap.put(it.key, conList = new LinkedList<String>())
            conList.addAll(it.value)
        }
    }
}

//
// CDR3 hypermutations
//

def cl2clCdr3Map = new HashMap<String, Map<String, Mutation>>()
def linkedByCdr3 = new HashSet<String>()

println "[${new Date()} $scriptName] Creating SHM links for CDR3"

clonotypeMap.clonotypes.each {
    def thisSeq = it.cdr3nt
    char[] thisNts = thisSeq.toCharArray()
    def fromClonotypes = clonotypeMap.getByCdr3(thisSeq)

    for (int i = 0; i < thisNts.length; i++) {
        char fromNt = thisNts[i]
        // hash-based 1-loop single-mm search
        Util.NTS.each { char toNt ->
            if (toNt != fromNt) {
                thisNts[i] = toNt

                def otherSeq = new String(thisNts)
                def toClonotypes = clonotypeMap.getByCdr3(otherSeq)

                if (toClonotypes && !linkedByCdr3.contains(thisSeq + "\t" + otherSeq)) {
                    // Make sure there are no duplicates
                    linkedByCdr3.add(otherSeq + "\t" + thisSeq)

                    int codonStart = i - i % 3

                    if (thisSeq.length() >= codonStart + 3) {
                        String fromCodon = thisSeq.substring(codonStart, codonStart + 3),
                               toCodon = otherSeq.substring(codonStart, codonStart + 3)

                        char fromAA = Util.codon2aa(fromCodon), toAA = Util.codon2aa(toCodon)

                        // Connect clonotypes with highest SHM overlap
                        Clonotype fromClonotype, toClonotype

                        (fromClonotype, toClonotype) = [fromClonotypes, toClonotypes].combinations().max {
                            Clonotype cl1 = it[0], cl2 = it[1]

                            cl1.shms.size() > cl2.shms.size() ?
                                    cl1.shms.intersect(cl2.shms).size() :
                                    cl2.shms.intersect(cl1.shms).size()
                        }

                        // Try to deduce SHM direction
                        int fromMutationSum = 0, toMutationSum = 0
                        if (deduceDirection) {
                            fromMutationSum = fromClonotypes.sum {
                                inDegreeMap[fromClonotype.key] ?: 0
                            } + fromClonotype.shms.size()
                            toMutationSum = toClonotypes.sum {
                                inDegreeMap[toClonotype.key] ?: 0
                            } + toClonotype.shms.size()
                        }

                        boolean flipped = false
                        if (toMutationSum < fromMutationSum) {
                            // flip
                            (fromNt, toNt) = [toNt, fromNt]
                            (fromAA, toAA) = [toAA, fromAA]
                            (fromClonotype, toClonotype) = [toClonotype, fromClonotype]
                            flipped = true
                        }

                        boolean directed = toMutationSum == fromMutationSum
                        def mutation = Mutation.cdr3Mutation(i, fromNt, toNt,
                                (int) (i / 3), fromAA, toAA, directed,
                                toClonotype)

                        if (!directed)
                            mutation.altFreq = (toClonotype.freq + fromClonotype.freq) / 2.0

                        def targetMap = cl2clCdr3Map[fromClonotype.key]

                        if (targetMap == null)
                            cl2clCdr3Map.put(fromClonotype.key,
                                    targetMap = new HashMap<String, Mutation>())

                        targetMap.put(toClonotype.key, mutation)

                        if (flipped)
                            fromNt = toNt // restore
                    }
                }
            }
            thisNts[i] = fromNt
        }
    }
}

//
// Write output
//

println "[${new Date()} $scriptName] Writing output"

def of = new File(outputPrefix).absoluteFile
if (of.parentFile != null)
    of.parentFile.mkdirs()

// Mutation lists

def alleles = clonotypeMap.clonotypes.collect { it.alleles }.flatten(),
    shms = clonotypeMap.clonotypes.collect { it.shms }.flatten(),
    shmsEmerged = [
            cl2clShmMap.values().collect { it.values() },
            cl2clCdr3Map.values().collect { it.values() }
    ].flatten()

// RS table

def allelesRsTable = new RSTable(true, vSegmentTable)
allelesRsTable.addAll(alleles)

def shmRsTable = new RSTable(true, vSegmentTable)
shmRsTable.addAll(shms)

def emergedShmRsTable = new RSTable(true, vSegmentTable)
emergedShmRsTable.addAll(shmsEmerged)

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

// Mutation motifs
def allelesPwm = new MotifPwm(), shmPwm = new MotifPwm(), emergedShmPwm = new MotifPwm()
allelesPwm.addAll(alleles)
shmPwm.addAll(shms)
emergedShmPwm.addAll(shmsEmerged)

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
    pw.println("key\t" + Clonotype.NODE_HEADER)
    clonotypeMap.clonotypes.each {
        pw.println("$it.key\t$it")
    }
}

new File(outputPrefix + ".edges.txt").withPrintWriter { pw ->
    pw.println("key\tweight\t" + Mutation.EDGE_HEADER)
    cl2clShmMap.each { from ->
        def fromKey = from.key
        from.value.each { to ->
            def toKey = to.key

            double weight = 1.0 / edgeWeightMap[toKey + "\t" + fromKey]

            to.value.each { shm ->
                pw.println(fromKey + " (" + shm.key + ") " + toKey + "\t" + weight + "\t" + shm)
            }
        }
    }
    cl2clCdr3Map.each { from ->
        from.value.each { to ->
            pw.println(from.key + " (CDR3) " + to.key + "\t" + 1.0 + "\t" + to.value)
        }
    }
}

new File(outputPrefix + ".net.txt").withPrintWriter { pw ->
    pw.println("source\tshm\ttarget")
    cl2clShmMap.each { from ->
        def fromKey = from.key
        from.value.each { to ->
            def toKey = to.key
            to.value.each { shm ->
                pw.println("$fromKey\t$shm.key\t$toKey")
            }
        }
    }
    cl2clCdr3Map.each {
        pw.println(it.key.replace(" (", "\t").replace(") ", "\t"))
    }
}

println "[${new Date()} $scriptName] Finished"