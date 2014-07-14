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

package com.antigenomics.vdjtools.intersection

import com.antigenomics.vdjtools.ClonotypeUtil
import com.antigenomics.vdjtools.Software
import com.antigenomics.vdjtools.spectratype.Spectratype
import org.apache.commons.io.FilenameUtils

def cli = new CliBuilder(usage: "IntersectSummary [options] sample1 sample2 output")
cli.h("display help message")
cli.S(longOpt: "software", argName: "string", required: true, args: 1,
        "Software used to process RepSeq data. Currently supported: ${Software.values().join(", ")}")

def opt = cli.parse(args)

if (opt == null)
    System.exit(-1)

if (opt.h || opt.arguments().size() < 3) {
    cli.usage()
    System.exit(-1)
}

def software = Software.byName(opt.S),
    sample1FileName = opt.arguments()[0], sample2FileName = opt.arguments()[1],
    outputFileName = opt.arguments()[2]

def scriptName = getClass().canonicalName.split("\\.")[-1]

//
// Load and intersect samples
//

println "[${new Date()} $scriptName] Reading samples"

def sample1 = ClonotypeUtil.loadClonotypes(sample1FileName, software),
    sample2 = ClonotypeUtil.loadClonotypes(sample2FileName, software)

println "[${new Date()} $scriptName] Intersecting"

def intersection = new PairedIntersection(sample1, sample2, IntersectionType.NucleotideV)

def intersectionResult = intersection.intersect(true)

//
// Write spectratypes for intersection
//

println "[${new Date()} $scriptName] Computing spectratypes"

def spec12 = new Spectratype(intersectionResult.clonotypes12),
    spec21 = new Spectratype(intersectionResult.clonotypes21)

println "[${new Date()} $scriptName] Writing output"

new File(outputFileName).withPrintWriter { pw ->

    // V segment, CDR3aa length, freq1, freq2, mean freq, correlation within set

    def header = "v\tlen\tn\t" +
            FilenameUtils.getBaseName(sample1FileName) + "_f\t" +
            FilenameUtils.getBaseName(sample2FileName) + "_f\t" +
            "F\tR"

    pw.println(header)

    double nonOverlappingFreq1 = intersectionResult.freq1 - intersectionResult.freq12,
           nonOverlappingFreq2 = intersectionResult.freq2 - intersectionResult.freq21

    // Overlapped and non-overlapped clonotypes

    def output = [
            ["Non-overlapping\t", intersectionResult.uniq1 + intersectionResult.uniq2 - intersectionResult.uniq12,
             nonOverlappingFreq1, nonOverlappingFreq2,
             Math.sqrt(nonOverlappingFreq1 * nonOverlappingFreq2),
             Double.NaN],

            ["Overlapping\t", intersectionResult.uniq12,
             intersectionResult.freq12, intersectionResult.freq21,
             intersectionResult.meanFrequency,
             intersectionResult.correlation]
    ]

    // Same for each spectratype peak separately

    output.addAll(
            (0..<spec12.numberOfPeaks).collect { int i ->
                def peak1 = spec12.sortedPeaks[i], peak2 = spec21.sortedPeaks[i]
                [peak1.signature, peak1.clones,
                 peak1.freq, peak2.freq,
                 Math.sqrt(peak1.freq * peak2.freq),
                 ClonotypeUtil.correlation(peak1.clonotypes, peak2.clonotypes)]
            }
    )

    output.sort { -it[3] }.each {
        pw.println(it.join("\t"))
    }

    /*pw12.println(spec12)
    pw12.println("Non-overlapping\t\t" +
            (intersectionResult.uniq1 - intersectionResult.uniq12) + "\t" +
            (intersectionResult.freq1 - intersectionResult.freq12))

    pw21.println(spec21)
    pw21.println("Non-overlapping\t\t" +
            (intersectionResult.uniq2 - intersectionResult.uniq12) + "\t" +
            (intersectionResult.freq2 - intersectionResult.freq21)) */
}

