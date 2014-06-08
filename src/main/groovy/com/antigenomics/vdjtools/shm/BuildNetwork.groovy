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

package com.antigenomics.vdjtools.shm

import com.antigenomics.vdjtools.Util

def cli = new CliBuilder(usage: 'BuildNetwork [options] input_level0 input_level1 input_level2 output_dir/')
cli.h('display help message')
def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 4) {
    cli.usage()
    System.exit(-1)
}

def nodeDataMap = new HashMap<String, String>(), edgeDataMap = new HashMap<String, String>()

String inputFileName1 = opt.arguments()[0], inputFileName2 = opt.arguments()[1], inputFileName3 = opt.arguments()[2],
       outputDir = opt.arguments()[3]

[inputFileName1, inputFileName2, inputFileName3].eachWithIndex { it, ind ->
    if (!new File(it).exists()) {
        println "[ERROR] Corresponding file ($it) for input level $ind does not exist."
    }
}

def NODE_HEADER = "key\tdisplay_name\tcount\t" +
        "cdr1nt\tcdr2nt\tcdr3nt\t" +
        "cdr1aa\tcdr2aa\tcdr3aa\t" +
        "v_segment\td_segment\tj_segment\t" +
        "cdr1q\tcdr2q\tcdr3q\t" +
        "mutations\t" +
        "inFrame\tnoStop\tcomplete"
// todo implement inframe, etc in IgBlastWrapper

def EDGE_HEADER = "key\tlevel\tcount\tnt_shm\taa_shm\tcdr_id\tsilent"

def l0Clonotypes = new HashSet<String>()

// 0    1   2   3                                               4   5   6               7           8           9
//17	.	.	TGTGTGAGACATAAACCTATGGTCCAGGGCGGCGTCGACGTCTGG	.	.	CVRHKPMVQGGVDVW	IGHV4-39*01	IGHD5-5*01	IGHJ6*01

def l0Key = { List<String> splitLine ->
    [splitLine[7..9], splitLine[3]].flatten().join("_")
}
def l1Key = { List<String> splitLine ->
    [splitLine[7..9], splitLine[1..3]].flatten().join("_")
}
def l2Key = { List<String> splitLine ->
    [splitLine[7..9], splitLine[1..3], splitLine[-1]].flatten().join("_")
}

new File(inputFileName1).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def nodeKey = l0Key(splitLine)
        nodeDataMap.put(nodeKey, splitLine[3] + "\t" + splitLine.join("\t"))
        l0Clonotypes.add(nodeKey)
    }
}

def level01Map = new HashMap<String, List>()

new File(inputFileName2).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def nodeKey = l1Key(splitLine),
            upperLevelKey = l0Key(splitLine)
        nodeDataMap.put(nodeKey, splitLine[1..3].join("-") + "\t" + splitLine.join("\t"))

        def otherLevel1Nodes = level01Map[upperLevelKey]
        if (!otherLevel1Nodes)
            level01Map.put(upperLevelKey, otherLevel1Nodes = new ArrayList())
        otherLevel1Nodes.add(nodeKey)
    }
}

def level12Map = new HashMap<String, List>()

new File(inputFileName3).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        if (splitLine[-1] != ".") {
            def nodeKey = l2Key(splitLine), upperLevelKey = l1Key(splitLine)
            nodeDataMap.put(nodeKey, "\t" + splitLine.join("\t"))

            def otherLevel2Nodes = level12Map[upperLevelKey]
            if (!otherLevel2Nodes)
                level12Map.put(upperLevelKey, otherLevel2Nodes = new ArrayList())
            otherLevel2Nodes.add([nodeKey, new HashSet(splitLine[-1].split("\\|").collect())])
        }
    }
}

def isSilent = { String shmString ->
    def fromToAA = shmString.split(",")[2].split(":")[1].split(">")
    fromToAA[0] == fromToAA[1]
}

// Report unique hypermutations only
level12Map.each {
    def upperLevelKey = it.key, l2family = it.value
    if (l2family.size() > 1) {  // if degree = 1 not worth reporting
        for (int i = 0; i < l2family.size(); i++) {
            def nodeId = l2family[i][0], shms = l2family[i][1]
            def shmIntersection = new HashSet()

            // Scan other L2 clonotypes with the same parent
            for (int j = i + 1; i < l2family.size(); i++) {
                shmIntersection.addAll(l2family[j][1].intersect(shms))
                // Remove intersection iteratively
                l2family[j][1].removeAll(shmIntersection)
            }

            // Remove all SHMs that intersect with otehr L2 clonotypes
            shms.removeAll(shmIntersection)

            // give all unique hypermutations their edges
            shms.each {
                edgeDataMap.put("$upperLevelKey (pp) $nodeId", "L2\t" +
                        it.replace(",", "\t") +
                        (it.contains("CDR") ? "" : "\t.") + "\t" +
                        isSilent(it))
            }
        }
    }
}

// Further decrease redundancy for L1 clonotypes
level01Map.each {
    def upperLevelKey = it.key, l0Degree = it.value.size()
    if (l0Degree > 1) // report all if other L1 clonotypes with the same parent are present
        it.value.each { nodeKey ->
            edgeDataMap.put("$upperLevelKey (pp) $nodeKey", "L1")
        }
    else
        it.value.each { nodeKey ->
            if (level12Map.containsKey(nodeKey) &&  // check for no hypermutations case
                    level12Map[nodeKey].size() > 1) // report if there are more than on L2 associated
                edgeDataMap.put("$upperLevelKey (pp) $nodeKey", "L1")
        }
}

// Hypermutations for CDR3 regions
def noShmL0Nodes = new LinkedList()
l0Clonotypes.each { thisNodeKey ->
    def splitKey = thisNodeKey.split("_")
    def vdjBase = splitKey[0..2].join("_"), cdr3Seq = splitKey[3]
    def chars = cdr3Seq.toCharArray()
    def oldNt
    boolean anyShms = false
    for (int i = 0; i < chars.length; i++) {
        oldNt = chars[i]
        // hash-based 1-loop single-mm search
        Util.NTS.each { char newNt ->
            if (newNt != oldNt) {
                chars[i] = newNt
                def otherSeq = new String(chars), otherNodeKey = vdjBase + "_" + otherSeq
                if (l0Clonotypes.contains(otherNodeKey) &&  // match exists
                        !edgeDataMap.containsKey("$otherNodeKey (pp) $thisNodeKey".toString()) // no duplicate edges
                ) {
                    int codonPos = i / 3
                    if (cdr3Seq.length() >= codonPos + 3) {
                        String thisCodon = cdr3Seq.substring(codonPos, codonPos + 3),
                               otherCodon = otherSeq.substring(codonPos, codonPos + 3)

                        def thisAA = Util.codon2aa(thisCodon), otherAA = Util.codon2aa(otherCodon),
                            silent = thisAA == otherAA

                        edgeDataMap.put("$thisNodeKey (pp) $otherNodeKey".toString(),
                                // todo: + cdr3 start
                                "L0\t.\t$i:$oldNt>$newNt\t$codonPos:$thisAA>$otherAA\tCDR3\t$silent")
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