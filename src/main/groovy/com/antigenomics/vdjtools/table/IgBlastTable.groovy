package com.antigenomics.vdjtools.table

import com.antigenomics.vdjtools.Util
import com.antigenomics.vdjtools.igblast.ClonotypeParsedData

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
def cli = new CliBuilder(usage: "IgBlastTable [options] igblast_output_level2 output_dir/")
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
       outputDir = opt.arguments()[1]

// Load regions
//IGHV1-18*01	1	75	76	99	100	150	151	174	175	288	VH	0
final Map<String, List<Integer>> regionSizesBySegment = new HashMap<>()

Util.resourceStreamReader("regions.${species}.txt").splitEachLine("\t") { List<String> splitLine ->
    regionSizesBySegment.put(splitLine[0],
            [splitLine[2].toInteger() - splitLine[1].toInteger() + 1,
             splitLine[4].toInteger() - splitLine[3].toInteger() + 1,
             splitLine[6].toInteger() - splitLine[5].toInteger() + 1,
             splitLine[8].toInteger() - splitLine[7].toInteger() + 1,
             splitLine[10].toInteger() - splitLine[9].toInteger() + 1])
}

println "[${new Date()} $scriptName] Looking at variable segments.."

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
        (1..6).collect { (0..<N_REGIONS).collect { Util.regionId2Name(it) }.join("\t") }.join('\t')

final Map<String, double[][][]> countersBySegmentByRegion = new HashMap<>(),
                                countersBySegmentByRegionPerBp = new HashMap<>()
//final Map<String, int[]> countersBySegment = new HashMap<>()

def allRegionCounters = new double[3][2][N_REGIONS], allRegionCountersPerBp = new double[3][2][N_REGIONS]
countersBySegmentByRegion.put("All", allRegionCounters)
countersBySegmentByRegionPerBp.put("All", allRegionCountersPerBp)
//countersBySegment.put("All", new int[2])

new File(inputFileNameL2).eachLine { line ->
    if (!line.startsWith("#")) {
        def clonotype = new ClonotypeParsedData(line)
        //def vSegment = splitLine[V_COL]
        //def vSegmentCounters = countersBySegment[vSegment]
        //if (vSegmentCounters == null)
        //    countersBySegment.put(vSegment, vSegmentCounters = new int[2])

        def regionSizes = regionSizesBySegment[clonotype.v]

        if (regionSizes == null) {
            println "[WARNING] No region info for $clonotype.v - skipping"
        } else {
            def vSegmentRegionCounters = countersBySegmentByRegion[clonotype.v],
                vSegmentRegionCountersPerBp = countersBySegmentByRegionPerBp[clonotype.v]

            if (vSegmentRegionCounters == null) {
                countersBySegmentByRegion.put(clonotype.v, vSegmentRegionCounters = new double[3][2][N_REGIONS])
                countersBySegmentByRegionPerBp.put(clonotype.v, vSegmentRegionCountersPerBp = new double[3][2][N_REGIONS])
            }

            clonotype.mutations.each { shmEntry ->
                int silentId = shmEntry.isSilent ? 1 : 0, regionId = Util.regionName2Id(shmEntry.region)

                double regionSize = regionSizes[regionId]

                def counters = [clonotype.count, clonotype.freq, 1],
                    perBpCounters = counters.collect { it / (double) regionSize }

                (0..2).each { int counterId ->
                    vSegmentRegionCounters[counterId][silentId][regionId] =
                            vSegmentRegionCounters[counterId][silentId][regionId] + counters[counterId]
                    allRegionCounters[counterId][silentId][regionId] =
                            allRegionCounters[counterId][silentId][regionId] + counters[counterId]

                    vSegmentRegionCountersPerBp[counterId][silentId][regionId] =
                            vSegmentRegionCountersPerBp[counterId][silentId][regionId] + perBpCounters[counterId]
                    allRegionCountersPerBp[counterId][silentId][regionId] =
                            allRegionCountersPerBp[counterId][silentId][regionId] + perBpCounters[counterId]
                }
            }
        }
    }
}

new File(outputDir).mkdirs()

new File(outputDir + "/mut_table.raw.txt").withPrintWriter { pw ->
    pw.println(HEADER)
    countersBySegmentByRegion.each {
        pw.println("$it.key\t${it.value.collect().flatten().join("\t")}")
    }
}

new File(outputDir + "/mut_table.perbp.txt").withPrintWriter { pw ->
    pw.println(HEADER)
    countersBySegmentByRegionPerBp.each {
        pw.println("$it.key\t${it.value.collect().flatten().join("\t")}")
    }
}