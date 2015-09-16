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


package com.antigenomics.vdjtools.util

import com.antigenomics.vdjtools.basic.SegmentUsage
import com.antigenomics.vdjtools.overlap.Overlap
import com.antigenomics.vdjtools.overlap.OverlapEvaluator
import com.antigenomics.vdjtools.sample.Sample

import java.nio.file.FileSystems
import java.nio.file.Path

import static java.io.File.separator

/**
 * Class that contains commonly used static functions for running scripts and I/O
 */
public class ExecUtil {
    // todo: use this everywhere
    public static final String MY_PATH = new File(ExecUtil.class.protectionDomain.codeSource.location.path).parent
    public static final int THREADS = Runtime.runtime.availableProcessors()

    /**
     * Gets the memory footprint of Java Runtime Environment
     * @return a string containing memory usage summary
     */
    public static String memoryFootprint() {
        final factor = 1024 * 1024 * 1024

        int maxMemory = Runtime.runtime.maxMemory() / factor,
            allocatedMemory = Runtime.runtime.totalMemory() / factor,
            freeMemory = Runtime.runtime.freeMemory() / factor

        "Memory usage: $allocatedMemory of ${maxMemory + freeMemory} GB"
    }

    public static void report(Object me, String message) {
        report(me, message, true)
    }

    public static void report(Object me, String message, boolean verbose) {
        if (verbose) {
            def scriptName = me.class.canonicalName.split("\\.")[-1]
            println "[${new Date()} $scriptName] $message"
        }
    }

    /**
     * Runs a specified Groovy script 
     * @param script groovy script object
     * @param args list of arguments
     * @return script
     */
    public static Object run(Script script, List<String> args) {
        // perform cleanup
        args = args.collect { it.trim() }
        args.removeAll { it.length() == 0 }
        //def argArray = args.split(" ").
        //        findAll { it != " " && it != "" }.
        //        collect { it.replaceAll("//+", "/").toString() }
        //println "Executing ${script.class.canonicalName} ${argArray.join(" ")}"
        //script.binding.setVariable("args", argArray)
        println "Executing ${script.class.canonicalName} ${args.collect().join(" ")}"
        script.binding.setVariable("args", args as String[])
        script.run()
    }

    /**
     * Routines for I/O logistics
     */

    /**
     * Makes sure that a given path exists
     * @param path
     */
    public static void ensureDir(String path) {
        if (path.endsWith(separator))
            new File(path).mkdirs()
        else
            new File(path).absoluteFile.parentFile.mkdirs()
    }

    /**
     * Creates a temporary directory
     * @param pathPrefix
     * @return
     */
    public static String createTempDir(String pathPrefix) {
        def tmpPath = getAbsolutePath(pathPrefix).toString() + "-vdjtools-" + UUID.randomUUID().toString()

        def tmpFolderFile = new File(tmpPath)
        tmpFolderFile.mkdirs()

        tmpFolderFile.deleteOnExit()

        tmpPath
    }

    /**
     * Internal
     */
    private static Path getAbsolutePath(String path) {
        getPath(path).toAbsolutePath()
    }

    /**
     * Internal
     */
    private static Path getPath(String path) {
        FileSystems.default.getPath(path)
    }

    /**
     * Provides a relative path to a sample given metadata path as reference.
     * Mainly used to store sample path in newly created metadata files
     * @param metadataPath path to a metadata file
     * @param samplePath path to a sample file
     * @return
     */
    public static String relativeSamplePath(String metadataPath, String samplePath) {
        getAbsolutePath(metadataPath).parent.relativize(getAbsolutePath(samplePath)).toString()
    }

    /**
     * Gets the absolute path to a sample file using metadata path as reference
     * @param metadataPath path to a metadata file
     * @param samplePath path to a sample file as recorded in metadata file
     * @return
     */
    public static String absoluteSamplePath(String metadataPath, String samplePath) {
        if (getPath(samplePath).absolute) {
            return samplePath
        }
        getAbsolutePath(metadataPath).parent.resolve(samplePath).normalize().toAbsolutePath().toString()
    }

    /**
     * Gets the output path according to VDJtools convention, with file name parts joined by '.'
     * @param outputPrefix output prefix, either directory name or directory name + prefix
     * @param outputSuffix one or more output suffices
     * @return output/prefix.suffix1.suffix2.txt or outputPrefix/suffix1.suffix2.txt
     */
    public static String formOutputPath(String outputPrefix,
                                        String... outputSuffix) {
        if (outputSuffix.any { it.contains(separator) })
            throw new IOException("Output suffices should not contain path separator")

        if (outputPrefix == ".")
            outputPrefix += separator

        boolean dir = outputPrefix.endsWith(separator)
        if (new File(outputPrefix).isDirectory() && !dir) {
            dir = true
            outputPrefix += separator
        }

        ensureDir(outputPrefix)

        def s = outputSuffix.join(".")

        (dir ? (outputPrefix + s) : (outputPrefix + "." + s)) +
                ((s.endsWith(".txt") || s.endsWith(".pdf") || s.endsWith(".png")) ? "" : ".txt")
    }

    public static String toPlotPath(String outputPath, String plotType) {
        outputPath[-4] == "." ? (outputPath[0..-5] + "." + plotType) : (outputPath + "." + plotType)
    }

    /**
     * Gets the output path for sample output according to VDJtools convention
     * @param outputPrefix output prefix, either directory name or directory name + prefix
     * @param sample sample object
     * @return
     */
    public static String formOutputPath(String outputPrefix, Sample sample) {
        formOutputPath(outputPrefix, sample.sampleMetadata.sampleId)
    }

    /**
     * Gets the metadata file output path according to VDJtools convention
     * @param outputPrefix output prefix, either directory name or directory name + prefix
     * @return
     */
    public static String formMetadataPath(String outputPrefix) {
        if (!new File(outputPrefix).isDirectory()) { // leave only directory in output prefix
            outputPrefix = getPath(outputPrefix).parent.toString()
        }

        formOutputPath(outputPrefix, "metadata")
    }

    public static void quiet() {
        SegmentUsage.VERBOSE = false
        OverlapEvaluator.VERBOSE = false
        Overlap.VERBOSE = false
    }
}
