package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.Util

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

def cli = new CliBuilder(usage: "IgBlastRsTable input output")
cli._(longOpt: 'no-norm', 'no normalization is performed (default is by length by v segment freq)')
def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() < 2) {
    cli.usage()
    System.exit(-1)
}

def scriptName = getClass().canonicalName.split("\\.")[-1]
boolean noNorm = opt.'no-norm'
String inputFileName = opt.arguments()[0],
       outputFileName = opt.arguments()[1]

int N_REGIONS = Util.N_REGIONS
final String HEADER = "Counter\t" +
        (0..<N_REGIONS).collect { "Frequency" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Frequency" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Clonotypes" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Clonotypes" }.join('\t') +
        "\nType\t" +
        (0..<N_REGIONS).collect { "Replacement" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Silent" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Replacement" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Silent" }.join('\t') +
        "\nRegion\t" +
        (1..4).collect { (0..<N_REGIONS).collect { Util.regionId2Name(it) }.join("\t") }.join('\t')

final Map<String, double[][][]> summaryTableBySegment = new HashMap<>()
def allTable =  new double[2][2][N_REGIONS]

summaryTableBySegment.put("All", allTable)


// 0        1       2           3               4       5       6       7       8           9               10      11
//#count	freq	v_segment	display_name	nt_shm	aa_shm	region	silent	undirected	region_length	v_freq	v_clonotypes

new File(inputFileName).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def (count, freq, vSegment, region, silent, regionLength, vFreq, vCount) =
          [splitLine[0].toDouble(), splitLine[1].toDouble(), splitLine[2],
          splitLine[6], splitLine[7].toBoolean() ? 1 : 0, splitLine[9].toInteger(),
          splitLine[10].toDouble(), splitLine[11].toDouble()]



        def regionId = Util.regionName2Id(region)
        def table = summaryTableBySegment[vSegment]
        if (table == null)
        summaryTableBySegment.put(vSegment, table = new double[2][2][N_REGIONS])

        if(!noNorm) {
            count /= (vCount * regionLength)
            freq /= (vFreq * regionLength)
        }

        table[0][silent][regionId] += count
        table[1][silent][regionId] += freq
        allTable[0][silent][regionId] += count
        allTable[1][silent][regionId] += freq
    }
}

new File(outputFileName).withPrintWriter { pw ->
    pw.println(HEADER)
    summaryTableBySegment.each {
        pw.println("$it.key\t${it.value.collect().flatten().join("\t")}")
    }
}