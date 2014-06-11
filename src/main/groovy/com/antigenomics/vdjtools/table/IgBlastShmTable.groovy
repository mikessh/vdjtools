package com.antigenomics.vdjtools.table

import com.antigenomics.vdjtools.RegionRanges
import com.antigenomics.vdjtools.igblast.MutationParseData

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

def FREQ_THRESHOLD = "0.40", SPEC_THRESHOLD = "3"
def cli = new CliBuilder(usage: "IgBlastNet [options] igblast_output_level2 output_dir/")
cli.h("display help message")
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

def freqThreshold = (opt.'allele-freq' ?: FREQ_THRESHOLD).toDouble(),
    spectraThreshold = (opt.'allele-spectra' ?: SPEC_THRESHOLD).toInteger()

String inputFileNameL2 = opt.arguments()[0],
       outputDir = opt.arguments()[1]

def V_COL = 13, CDR_COL = 4..6, CDR3_COL = 6, SHM_COL = -1,
    CDR3_AA_COL = 9, EVENT_FREQ_COL = 3

final Map<String, double[][][]> countersBySegmentByRegion = new HashMap<>()
final Map<String, int[]> countersBySegment = new HashMap<>()

println "[${new Date()} $scriptName] Looking at variable segments.."

int N_REGIONS = 7
final String HEADER = "Counter\t" +
        (0..<N_REGIONS).collect { "Reads" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Reads" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Ratio" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Ratio" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Clones" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Clones" }.join('\t') +
        "\nType\t" +
        (0..<N_REGIONS).collect { "Replacement" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Silent" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Replacement" }.join('\t') + '\t' +
        (0..<N_REGIONS).collect { "Silent" }.join('\t') +
        "\nChain\t" +
        (0..3).collect { RegionRanges.HEADER }.join('\t')

new File(inputFileNameL2).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def vSegment = splitLine[V_COL]
        def vSegmentConters = countersBySegment[vSegment]
        if (vSegmentConters == null) {
            countersBySegment.put(vSegment, new int[2])
            countersBySegmentByRegion.put(vSegment, new double[3][2][N_REGIONS])
        }
    }
}

countersBySegmentByRegion.put("All", new double[3][2][5])
countersBySegment.put("All", new int[2])

new File(inputFileNameL2).splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def vSegment = splitLine[V_COL]
        def vSegmentConters = countersBySegment[vSegment]
        if (vSegmentConters == null) {
            countersBySegment.put(vSegment, vSegmentConters = new int[2])
            countersBySegmentByRegion.put(it[0], new double[3][2][5])
            countersBySegment.put(it[0], new int[2])
        }
        splitLine[SHM_COL].split("\\|").collect { new MutationParseData(it) }.each { shmEntry ->

        }
    }
}