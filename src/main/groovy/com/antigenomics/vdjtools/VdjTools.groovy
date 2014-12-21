/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 *
 * Last modified on 8.11.2014 by mikesh
 */


package com.antigenomics.vdjtools

import com.antigenomics.vdjtools.basic.*
import com.antigenomics.vdjtools.db.ScanDatabase
import com.antigenomics.vdjtools.diversity.CalcDiversityStats
import com.antigenomics.vdjtools.diversity.DownSample
import com.antigenomics.vdjtools.diversity.RarefactionPlot
import com.antigenomics.vdjtools.intersection.*
import com.antigenomics.vdjtools.pwm.ComputePwms
import com.antigenomics.vdjtools.util.ExecUtil
import com.antigenomics.vdjtools.util.RInstall

import java.util.jar.JarFile

def version = (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
        /Implementation-Version: (.+)/)[0][1]

def printHelp = {
    println "VDJtools V$version"
    println ""
    println "Run as \$java -jar vdjtools-${version}.jar ROUTINE_NAME arguments"
    println ""
    println "[Single-sample statistics]"
    println "CalcBasicStats"
    println "CalcSpectratype"
    println "CalcSegmentUsage"
    println "PlotFancySpectratype"
    println "PlotSpectratypeV"
    println "PlotFancyVJUsage"
    println "FilterNonFunctional"
    println ""
    println "[Sample diversity]"
    println "CalcDiversityStats"
    println "RarefactionPlot"
    println "DownSample"
    //println "BuildFrequencyTable"
    println ""
    println "[Cross-sample analysis]"
    println "IntersectPair"
    println "BatchIntersectPair"
    println "BatchIntersectPairPlot"
    println "IntersectSequential"
    println "PoolSamples"
    println ""
    println "[Sample annotation]"
    println "ScanDatabase"
    println "ComputePwms"
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
        case "FILTERNONFUNCTIONAL":
            return new FilterNonFunctional()
        case "DOWNSAMPLE":
            return new DownSample()
        case "INTERSECTPAIR":
            return new IntersectPair()
        case "BATCHINTERSECTPAIR":
            return new BatchIntersectPair()
        case "BATCHINTERSECTPAIRPLOT":
            return new BatchIntersectPairPlot()
        case "INTERSECTSEQUENTIAL":
            return new IntersectSequential()
        case "SCANDATABASE":
            return new ScanDatabase()
    //case "BUILDFREQUENCYTABLE":
    //    return new BuildFrequencyTable()
        case "POOLSAMPLES":
            return new PoolSamples()
        case "RAREFACTIONPLOT":
            return new RarefactionPlot()
        case "COMPUTEPWMS":
            return new ComputePwms()
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
        ExecUtil.run(script, args.length > 1 ? args[1..-1].join(" ") : "")
    } catch (Exception e) {
        println "[ERROR] $e.message, see _vdjtools_error.log for details"
        new File("_vdjtools_error.log").withWriterAppend { writer ->
            writer.println("[${new Date()}]")
            writer.println("[Script]")
            writer.println(args[0])
            writer.println("[CommandLine]")
            writer.println("executing vdjtools-${version}.jar ${args.join(" ")}")
            writer.println("[Message]")
            writer.println(e.message)
            writer.println("[StackTrace-Short]")
            writer.println(e.stackTrace.findAll { it.toString().contains("com.antigenomics.vdjtools") }.join("\n"))
            writer.println("[StackTrace-Full]")
            e.printStackTrace(new PrintWriter(writer))
        }
        System.exit(-1)
    }
}