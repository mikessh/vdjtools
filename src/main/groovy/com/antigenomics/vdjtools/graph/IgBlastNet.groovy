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

def FREQ_THRESHOLD = "0.40", SPEC_THRESHOLD = "3"
def cli = new CliBuilder(usage: "IgBlastNet [options] igblast_output_level0 igblast_output_level2 output_dir/")
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

String inputFileNameL0 = opt.arguments()[0],
       inputFileNameL2 = opt.arguments()[1],
       outputDir = opt.arguments()[2]

[(inputFileNameL0): 0, (inputFileNameL2): 2].each {
    if (!new File(it.key).exists()) {
        println "[ERROR] Corresponding file ($it.key) for input level $it.value does not exist."
    }
}

//
// Parsing
//

def V_COL = 13, CDR_COL = 4..6, CDR3_COL = 6, SHM_COL = -1,
    CDR3_AA_COL = 9, EVENT_FREQ_COL = 3

def l0Key = { List<String> splitLine ->
    [splitLine[V_COL], splitLine[CDR3_COL]].flatten().join("_")
}

def l2Key = { List<String> splitLine ->
    [splitLine[V_COL], splitLine[CDR_COL], splitLine[SHM_COL] == "." ? "." :
            splitLine[SHM_COL].split("\\|").collect {
                it.split(",")[1]
            }].flatten().join("_")
}

//
// Network storage
//

def nodeDataMap = new HashMap<String, String>(), edgeDataMap = new HashMap<String, String>()

//
// Load level 0 clonotypes, unique CDR3
//

println "[${new Date()} $scriptName] Loading L0 clonotypes"

def l0Clonotypes = new HashSet<String>()
new File(inputFileNameL0).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def nodeKey = l0Key(splitLine)
        nodeDataMap.put(nodeKey, "L0\t" + splitLine[V_COL][3] + splitLine[V_COL][2] +
                splitLine[V_COL][4..-4] + ":" + splitLine[CDR3_AA_COL] + "\t"
                + splitLine.join("\t"))
        l0Clonotypes.add(nodeKey)
    }
}

//
// Storage for SHM data
//

class ShmParsedData {
    final String ntString, aaString, region, countString
    final ntPos, aaPos
    final char fromAA, toAA, fromNT, toNT
    final boolean isSilent

    ShmParsedData(String shmString) {
        def splitString = shmString.split(",")
        countString = splitString[0].replace(":", "\t")

        ntString = splitString[1]
        def splitNTString = ntString.split("[:>]")
        ntPos = splitNTString[0].toInteger()
        fromNT = splitNTString[1]
        toNT = splitNTString[2]

        aaString = splitString[2]
        def splitAAString = aaString.split("[:>]")
        aaPos = splitAAString[0].toInteger()
        fromAA = splitAAString[1]
        toAA = splitAAString[2]
        isSilent = fromAA == toAA

        region = splitString[3]
    }

    String edgeString() {
        ["L2", isSilent ? "S" : (region + ":" + fromAA + ">" + toAA),
         countString, ntString, aaString, region, isSilent].join("\t")
    }


    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        ShmParsedData that = (ShmParsedData) o

        if (ntString != that.ntString) return false

        return true
    }

    int hashCode() {
        return ntString.hashCode()
    }
}

class ShmGlobalCounter {
    double freq = 0
    final Set<Integer> cdr3Len = new HashSet<>()
}

def shmFreqBySegmentMap = new HashMap<String, Map<ShmParsedData, ShmGlobalCounter>>()
def vSegmFreqMap = new HashMap<String, Double>()

//
// Load level 2 clonotypes, unique CDR3+V+mutations
//

println "[${new Date()} $scriptName] Loading L2 clonotypes"

def level02Map = new HashMap<String, List>()
new File(inputFileNameL2).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        if (splitLine[SHM_COL] != ".") {
            def nodeKey = l2Key(splitLine), upperLevelKey = l0Key(splitLine)
            nodeDataMap.put(nodeKey, "L2\t\t" + splitLine.join("\t"))

            def otherLevel2Nodes = level02Map[upperLevelKey]
            if (!otherLevel2Nodes)
                level02Map.put(upperLevelKey, otherLevel2Nodes = new ArrayList())

            def shms = new HashSet(splitLine[-1].split("\\|").collect { new ShmParsedData(it) })
            otherLevel2Nodes.add([nodeKey, shms])

            def freq = splitLine[EVENT_FREQ_COL].toDouble()
            def vSegment = splitLine[V_COL]

            vSegmFreqMap.put(vSegment, (vSegmFreqMap[vSegment] ?: 0) + freq)
            def shmFreqMap = shmFreqBySegmentMap[vSegment]
            if (shmFreqMap == null)
                shmFreqBySegmentMap.put(vSegment, shmFreqMap = new HashMap<ShmParsedData, ShmGlobalCounter>())

            shms.each { ShmParsedData shmData ->
                def shmCounter = shmFreqMap[shmData]
                if (shmCounter == null)
                    shmFreqMap.put(shmData, shmCounter = new ShmGlobalCounter())
                shmCounter.freq += freq
                shmCounter.cdr3Len.add(splitLine[CDR3_AA_COL].length())
            }
        }
    }
}

//
// Report hypermutations, check for redundancy
//

println "[${new Date()} $scriptName] Filtering mutations"

level02Map.each {
    def upperLevelKey = it.key, l2family = it.value
    if (l2family.size() > 1) {
        // if degree > 1, we just remove all redundant mutations
        for (int i = 0; i < l2family.size(); i++) {
            def nodeId = l2family[i][0]
            Set<ShmParsedData> shms = l2family[i][1]
            def shmIntersection = new HashSet()

            // Scan other L2 clonotypes with the same parent
            for (int j = i + 1; i < l2family.size(); i++) {
                Set otherShms = l2family[j][1]
                shmIntersection.addAll(otherShms.intersect(shms))
                // Remove intersection iteratively
                otherShms.removeAll(shmIntersection)
            }

            // Remove all SHMs that intersect with other L2 clonotypes
            shms.removeAll(shmIntersection)

            // Give unique hypermutations their edges, check for alleles
            shms.each {
                edgeDataMap.put("$upperLevelKey (pp) $nodeId", it.edgeString())
            }
        }
    } else {
        // in case of single 2nd level clonotype,
        // we just check if the mutation is germline using an empirical rule
        def vSegment = upperLevelKey.split("_")[0]
        def nodeId = l2family[0][0]
        Set<ShmParsedData> shms = l2family[0][1]
        def shmFreqMap = shmFreqBySegmentMap[vSegment], vFreq = vSegmFreqMap[vSegment]
        shms.each {
            def shmCounters = shmFreqMap[it]
            if (shmCounters.freq / vFreq >= freqThreshold && shmCounters.cdr3Len.size() >= spectraThreshold)
                edgeDataMap.put("$upperLevelKey (pp) $nodeId", it.edgeString())
        }
    }
}

//
// Hypermutations for CDR3 regions
//

println "[${new Date()} $scriptName] Looking for CDR3 mismatch pairs"

def noShmL0Nodes = new LinkedList()
l0Clonotypes.each { thisNodeKey ->
    def splitKey = thisNodeKey.split("_")
    def vSegment = splitKey[0], cdr3Seq = splitKey[1]
    def chars = cdr3Seq.toCharArray()
    def oldNt
    boolean anyShms = false
    for (int i = 0; i < chars.length; i++) {
        oldNt = chars[i]
        // hash-based 1-loop single-mm search
        Util.NTS.each { char newNt ->
            if (newNt != oldNt) {
                chars[i] = newNt
                def otherSeq = new String(chars), otherNodeKey = vSegment + "_" + otherSeq
                if (l0Clonotypes.contains(otherNodeKey) &&  // match exists
                        !edgeDataMap.containsKey("$otherNodeKey (pp) $thisNodeKey".toString()) // no duplicate edges
                ) {
                    int codonStart = i - i % 3
                    if (cdr3Seq.length() >= codonStart + 3) {
                        String thisCodon = cdr3Seq.substring(codonStart, codonStart + 3),
                               otherCodon = otherSeq.substring(codonStart, codonStart + 3)

                        def thisAA = Util.codon2aa(thisCodon), otherAA = Util.codon2aa(otherCodon),
                            silent = thisAA == otherAA

                        edgeDataMap.put("$thisNodeKey (pp) $otherNodeKey".toString(),
                                "L0\t${silent ? "S" : "CDR3:$thisAA<>$otherAA"}\t" +
                                        ".\t.\t.\t.\t$i:$oldNt<>$newNt\t${(int) (i / 3)}:$thisAA>$otherAA\tCDR3\t$silent")
                        anyShms = true
                    }
                }
            }
        }
        chars[i] = oldNt
    }
    if (!anyShms)
        noShmL0Nodes.add(thisNodeKey)
}

//
// Output
//

def NODE_HEADER = "key\tlevel\tdisplay_name\treads\tfreq_reads\tevents\tfreq_events\t" +
        "cdr1nt\tcdr2nt\tcdr3nt\t" +
        "cdr1aa\tcdr2aa\tcdr3aa\t" +
        "inFrame\tnoStop\tcomplete\t" +
        "v_segment\td_segment\tj_segment\t" +
        "cdr1q\tcdr2q\tcdr3q\t" +
        "mutations"

def EDGE_HEADER = "key\tlevel\tdisplay_name\t" +
        "count_reads\tfreq_reads\tcount_events\tfreq_events\t" +
        "nt_shm\taa_shm\tregion\tsilent"

println "[${new Date()} $scriptName] Writing output"

new File(outputDir).mkdirs()

new File("$outputDir/nodes.txt").withPrintWriter { pw ->
    pw.println(NODE_HEADER)
    nodeDataMap.each {
        pw.println("$it.key\t$it.value")
    }
}

new File("$outputDir/edges.txt").withPrintWriter { pw ->
    pw.println(EDGE_HEADER)
    edgeDataMap.each {
        pw.println("$it.key\t$it.value")
    }
}

new File("$outputDir/net.txt").withPrintWriter { pw ->
    pw.println("source\ttarget")
    noShmL0Nodes.each {
        pw.println("$it\t")
    }
    edgeDataMap.each {
        pw.println(it.key.replace(" (pp) ", "\t"))
    }
}

println "[${new Date()} $scriptName] Finished"