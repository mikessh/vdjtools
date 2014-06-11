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

package com.antigenomics.vdjtools.graph

import com.antigenomics.vdjtools.Util
import com.antigenomics.vdjtools.igblast.ClonotypeParsedData
import com.antigenomics.vdjtools.igblast.MutationParseData

def FREQ_THRESHOLD = "0.40", SPEC_THRESHOLD = "3"
def cli = new CliBuilder(usage: "IgBlastNet [options] igblast_output_level2 output_dir/")
cli.h("display help message")
cli._(longOpt: "allele-freq", argName: "[0, 1]", "Frequency threshold, used together with spectratype threshold. " +
        "Mutations with higher frequency are considered as allele candidates. [default=$FREQ_THRESHOLD]")
cli._(longOpt: "allele-spectra", argName: "int", "Spectratype threshold, used together with frequency threshold. " +
        "At least \$allele-spectra clonotypes with distinct CDR3 lengths must contain this mutation " +
        "for it to be considered as allele. [default=$SPEC_THRESHOLD]")
def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 3) {
    cli.usage()
    System.exit(-1)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]

def freqThreshold = (opt.'allele-freq' ?: FREQ_THRESHOLD).toDouble(),
    spectraThreshold = (opt.'allele-spectra' ?: SPEC_THRESHOLD).toInteger()

String inputFileNameL2 = opt.arguments()[1],
       outputDir = opt.arguments()[2]

if (!new File(inputFileNameL2).exists()) {
    println "[ERROR] Input file does not exist: $inputFileNameL2"
}

//
// Allele statistics
//

class MutCounter {
    double freq = 0
    final Set<Integer> cdr3Len = new HashSet<>()
}

def freqByV = new HashMap<String, Double>(),
    freqByMutByV = new HashMap<String, Map<MutationParseData, MutCounter>>()

//
// Read and parse clonotypes
//

def clonotypeMap = new HashMap<String, ClonotypeParsedData>(),
    byCdr3Map = new HashMap<String, List<ClonotypeParsedData>>()

println "[${new Date()} $scriptName] Loading L2 clonotypes"

new File(inputFileNameL2).eachLine { line ->
    if (!line.startsWith("#")) {
        def clonotype = new ClonotypeParsedData(line)
        clonotypeMap.put(clonotype.key, clonotype)

        freqByV.put(clonotype.v, (freqByV[clonotype.v] ?: 0) + clonotype.freq)

        def freqByMut = freqByMutByV[clonotype.v]

        if (freqByMut == null)
            freqByMutByV.put(clonotype.v, freqByMut = new HashMap<MutationParseData, MutCounter>())

        clonotype.mutations.each { mpd ->
            def mc = freqByMut[mpd]
            if (mc == null)
                freqByMut.put(mpd, mc = new MutCounter())
            mc.freq += clonotype.freq
            mc.cdr3Len.add(clonotype.cdr3nt.length())
        }

        def clonotypeList = byCdr3Map[clonotype.cdr3nt]
        if (clonotypeList == null)
            byCdr3Map.put(clonotype.cdr3nt, clonotypeList = new ArrayList<ClonotypeParsedData>())
        clonotypeList.add(clonotype)
    }
}

//
// Deduce alleles
//

println "[${new Date()} $scriptName] Deducing alleles"

clonotypeMap.values().each { clonotype ->
    def vFreq = freqByV[clonotype.v],
        freqByMut = freqByMutByV[clonotype.v]

    clonotype.mutations.each { mpd ->
        def mutCounter = freqByMut[mpd]
        if (mutCounter.freq / vFreq < freqThreshold ||
                mutCounter.cdr3Len.size() < spectraThreshold)
            clonotype.shms.add(mpd)
        else
            clonotype.alleles.add(mpd)
    }
}

//
// Iterate through clonotypes with same CDR3nt and build hypermutation links
//

def cl2clShmMap = new HashMap<String, Map<String, List<MutationParseData>>>()

def putShm = { String key1, String key2, MutationParseData shm ->
    def clShmMap = cl2clShmMap[key2]
    if (clShmMap == null)
        cl2clShmMap.put(key2,
                clShmMap = new HashMap<String, List<MutationParseData>>())
    def shmList = clShmMap[key1]
    if (shmList == null)
        clShmMap.put(key1, shmList = new LinkedList<MutationParseData>())
    shmList.add(shm)
}

class Checker {
    // Recursively scan graph
    static boolean checkNodes(HashMap<String, List<String>> graph, String clonotype1, String clonotype2) {
        def subNetwork = graph[clonotype1]
        subNetwork.each {
            if (it == clonotype2)
                return true
            else if (checkNodes(graph, it, clonotype2))
                return true
        }
        return false
    }
}

println "[${new Date()} $scriptName] Creating SHM links for FW1-FW3"

byCdr3Map.values().each { clonotypeList ->
    int maxLevel = 0

    clonotypeList.each { maxLevel = Math.max(maxLevel, it.shms.size()) }

    final def connectivityMap = new HashMap<String, List<String>>(),
              newConnectivityMap = new HashMap<String, List<String>>()

    def addClonotypePair = { String key1, String key2 ->
        def conList = newConnectivityMap[key2]
        if (conList == null)
            newConnectivityMap.put(key2, conList = new LinkedList<String>())
        conList.add(key1)
        conList = newConnectivityMap[key1]
        if (conList == null)
            newConnectivityMap.put(key1, conList = new LinkedList<String>())
        conList.add(key2)
    }

    def checkAndAddShm = { ClonotypeParsedData clonotype1, ClonotypeParsedData clonotype2, MutationParseData shm ->
        def key1 = clonotype1.key, key2 = clonotype2.key
        // check if not connected yet
        if (!Checker.checkNodes(connectivityMap, key1, key2)) {
            putShm(key1, key2, shm)

            // update connectivity map
            addClonotypePair(key1, key2)
        }
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
                if (Math.abs(shmDelta1) + Math.abs(shmDelta2) == level) {
                    shms1.each { shm ->
                        if (!shmIntersection.contains(shm))
                            checkAndAddShm(clonotype1, clonotype2, shm)
                    }
                    shms2.each { shm ->
                        if (!shmIntersection.contains(shm))
                            checkAndAddShm(clonotype2, clonotype1, shm)
                    }
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

//
// CDR3 hypermutations
//

def cl2clCdr3Map = new HashMap<String, String>()

println "[${new Date()} $scriptName] Creating SHM links for CDR3"

clonotypeMap.values().each {
    def thisSeq = it.cdr3nt
    def ntChars = thisSeq.toCharArray()

    for (int i = 0; i < ntChars.length; i++) {
        oldNt = ntChars[i]
        // hash-based 1-loop single-mm search
        Util.NTS.each { char newNt ->
            if (newNt != oldNt) {
                ntChars[i] = newNt
                def otherSeq = new String(ntChars)
                if (byCdr3Map.containsKey(otherSeq)) {
                    if (!cl2clCdr3Map.containsKey(otherSeq + " (CDR3) " + thisSeq)) {
                        int codonStart = i - i % 3
                        if (thisSeq.length() >= codonStart + 3) {
                            String thisCodon = thisSeq.substring(codonStart, codonStart + 3),
                                   otherCodon = otherSeq.substring(codonStart, codonStart + 3)

                            def thisAA = Util.codon2aa(thisCodon), otherAA = Util.codon2aa(otherCodon),
                                silent = thisAA == otherAA
                            cl2clCdr3Map.put(thisSeq + " (CDR3) " + otherSeq,
                                    // display name
                                    "${silent ? "S" : "CDR3:$thisAA<>$otherAA"}\t" +
                                            // tostring
                                            "$i:$oldNt<>$newNt\t${(int) (i / 3)}:$thisAA>$otherAA\tCDR3\t$silent")
                        }
                    }
                }
            }
            ntChars[i] = oldNt
        }
    }
}

//
// Write output
//

println "[${new Date()} $scriptName] Writing output"

new File(outputDir).mkdirs()

new File("$outputDir/nodes.txt").withPrintWriter { pw ->
    pw.println("key\t" + ClonotypeParsedData.NODE_HEADER)
    clonotypeMap.each {
        pw.println("$it.key\t$it.value")
    }
}

new File("$outputDir/edges.txt").withPrintWriter { pw ->
    pw.println("key\t" + MutationParseData.EDGE_HEADER)
    cl2clShmMap.each { from ->
        def fromKey = from.key
        from.value.each { to ->
            def toKey = to.key
            to.value.each { shm ->
                pw.println(fromKey + " (" + shm.key + ") " + toKey + "\t" + shm)
            }
        }
    }
    cl2clCdr3Map.each {
        pw.println("$it.key\t$it.value")
    }
}

new File("$outputDir/net.txt").withPrintWriter { pw ->
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