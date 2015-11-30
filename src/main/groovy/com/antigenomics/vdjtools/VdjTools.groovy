/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */


package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.basic.*
import com.antigenomics.vdjtools.compare.Enrichment
import com.antigenomics.vdjtools.diversity.CalcDiversityStats
import com.antigenomics.vdjtools.diversity.PlotQuantileStats
import com.antigenomics.vdjtools.diversity.RarefactionPlot
import com.antigenomics.vdjtools.operate.JoinSamples
import com.antigenomics.vdjtools.operate.PoolSamples
import com.antigenomics.vdjtools.overlap.*
import com.antigenomics.vdjtools.preprocess.*
import com.antigenomics.vdjtools.profile.CalcCdrAAProfile
import com.antigenomics.vdjtools.util.*

import java.util.jar.JarFile

def version = (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
        /Implementation-Version: (.+)/)[0][1]

def printHelp = {
    println "VDJtools V$version"
    println ""
    println "Run as \$java -jar vdjtools-${version}.jar ROUTINE_NAME arguments"
    println ""
    println "[Basic]"
    println "CalcBasicStats"
    println "CalcSpectratype"
    println "CalcSegmentUsage"
    println "PlotFancySpectratype"
    println "PlotSpectratypeV"
    println "PlotFancyVJUsage"
    println ""
    println "[Diversity]"
    println "CalcDiversityStats"
    println "RarefactionPlot"
    println "PlotQuantileStats"
    println ""
    println "[Overlap]"
    println "OverlapPair"
    println "CalcPairwiseDistances"
    println "ClusterSamples"
    println "TestClusters"
    println "TrackClonotypes"
    println ""
    println "[Preprocessing]"
    println "ApplySampleAsFilter"
    println "FilterNonFunctional"
    println "DownSample"
    println "Decontaminate"
    println "FilterBySegment"
    println "SelectTop"
    println ""
    println "[Operation]"
    println "PoolSamples"
    println "JoinSamples"
    println "Enrichment"
    println ""
    println "[Annotation]"
    println "(ScanDatabase) -> moved to VDJdb since 1.0.5"
    println "CalcCdrAAProfile"
    println ""
    println "[Util]"
    println "FilterMetadata"
    println "SplitMetadata"
    println "Convert"
    println "RInstall"
}

def getScript = { String scriptName ->
    switch (scriptName.toUpperCase()) {
        case "CALCBASICSTATS":
            return new CalcBasicStats()
        case "CALCSPECTRATYPE":
            return new CalcSpectratype()
        case "CALCSEGMENTUSAGE":
            return new CalcSegmentUsage()
        case "PLOTFANCYSPECTRATYPE":
            return new PlotFancySpectratype()
        case "PLOTSPECTRATYPEV":
            return new PlotSpectratypeV()
        case "PLOTFANCYVJUSAGE":
            return new PlotFancyVJUsage()

        case "CALCDIVERSITYSTATS":
            return new CalcDiversityStats()
        case "RAREFACTIONPLOT":
            return new RarefactionPlot()
        case "PLOTQUANTILESTATS":
            return new PlotQuantileStats()

        case "OVERLAPPAIR":
            return new OverlapPair()
        case "CALCPAIRWISEDISTANCES":
            return new CalcPairwiseDistances()
        case "CLUSTERSAMPLES":
            return new ClusterSamples()
        case "TESTCLUSTERS":
            return new TestClusters()
        case "TRACKCLONOTYPES":
            return new TrackClonotypes()

        case "APPLYSAMPLEASFILTER":
            return new ApplySampleAsFilter()
        case "DECONTAMINATE":
            return new Decontaminate()
        case "FILTERBYSEGMENT":
            return new FilterBySegment()
        case "FILTERNONFUNCTIONAL":
            return new FilterNonFunctional()
        case "DOWNSAMPLE":
            return new DownSample()
        case "SELECTTOP":
            return new SelectTop()
        case "CORRECT":
            return new Correct()

        case "POOLSAMPLES":
            return new PoolSamples()
        case "JOINSAMPLES":
            return new JoinSamples()
        case "ENRICHMENT":
            return new Enrichment()

        case "CALCCDRAAPROFILE":
            return new CalcCdrAAProfile()
        case "SCANDATABASE":
            println "Moved to VDJdb since 1.0.5, see docs"
            System.exit(0)
            break

        case "FILTERMETADATA":
            return new FilterMetadata()
        case "SPLITMETADATA":
            return new SplitMetadata()
        case "CONVERT":
            return new Convert()
        case "RINSTALL":
            return new RInstall()
        case "-H":
        case "H":
        case "-HELP":
        case "HELP":
        case "":
            printHelp()
            println ""
            System.exit(0)
            break

        default:
            printHelp()
            println ""
            println "Unknown routine name $scriptName"
            System.exit(0)
    }
}

if (args.length == 0)
    printHelp()
else {
    def script = getScript(args[0])
    try {
        ExecUtil.run(script, args.length > 1 ? args[1..-1] : [""])
    } catch (Exception e) {
        println "[ERROR] ${e.toString()}, see _vdjtools_error.log for details"
        new File("_vdjtools_error.log").withWriterAppend { writer ->
            writer.println("[${new Date()} BEGIN]")
            writer.println("[Script]")
            writer.println(args[0])
            writer.println("[CommandLine]")
            writer.println("executing vdjtools-${version}.jar ${args.join(" ")}")
            writer.println("[Message]")
            writer.println(e.toString())
            writer.println("[StackTrace-Short]")
            writer.println(e.stackTrace.findAll { it.toString().contains("com.antigenomics.vdjtools") }.join("\n"))
            writer.println("[StackTrace-Full]")
            e.printStackTrace(new PrintWriter(writer))
            writer.println("[END]")
        }
        System.exit(1)
    }
}