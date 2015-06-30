/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
import com.antigenomics.vdjtools.util.Convert
import com.antigenomics.vdjtools.util.ExecUtil
import com.antigenomics.vdjtools.util.RInstall
import com.antigenomics.vdjtools.util.ScanDatabase

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
    println ""
    println "[Operation]"
    println "PoolSamples"
    println "JoinSamples"
    println "Enrichment"
    println ""
    println "[Annotation]"
    println "ScanDatabase"
    println "CalcCdrAAProfile"
    println ""
    println "[Util]"
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

        case "POOLSAMPLES":
            return new PoolSamples()
        case "JOINSAMPLES":
            return new JoinSamples()
        case "ENRICHMENT":
            return new Enrichment()

        case "SCANDATABASE":
            return new ScanDatabase()
        case "CALCCDRAAPROFILE":
            return new CalcCdrAAProfile()

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
            System.exit(-1)
            break

        default:
            printHelp()
            println ""
            println "Unknown routine name $scriptName"
            System.exit(-1)
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
        System.exit(-1)
    }
}