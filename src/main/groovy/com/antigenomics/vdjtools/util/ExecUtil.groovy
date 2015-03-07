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


package com.antigenomics.vdjtools.util

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
                ((s.endsWith(".txt") || s.endsWith(".pdf")) ? "" : ".txt")
    }

    public static String toPlotPath(String outputPath) {
        outputPath[-4] == "." ? (outputPath[0..-5] + ".pdf") : (outputPath + ".pdf")
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
}
