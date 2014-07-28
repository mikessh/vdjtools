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

package com.antigenomics.vdjtools.util

import com.antigenomics.vdjtools.sample.metadata.MetadataEntry

class RUtil {
    static String asNumeric(MetadataEntry metadataEntry) {
        def value = metadataEntry.asNumeric()
        value.isNaN() ? "NA" : value.toString()
    }

    static void execute(String scriptName, String... params) {
        // Create a temp file to store the script
        def scriptRes = CommonUtil.resourceStreamReader("rscripts/$scriptName")
        scriptName = UUID.randomUUID().toString() + "_" + scriptName

        def scriptFile = new File(scriptName)

        scriptFile.withPrintWriter { pw ->
            scriptRes.readLines().each {
                pw.println(it)
            }
        }

        scriptFile.deleteOnExit()

        // Run script
        def cmd = ["Rscript", scriptName, params]
        println "[RUtil] Executing ${cmd.flatten().join(" ")}"
        def proc = cmd.flatten().execute()

        proc.in.eachLine {
            println(it)
        }

        proc.out.close()
        proc.waitFor()

        if (proc.exitValue()) {
            println "[ERROR] ${proc.getErrorStream()}"
        }
    }

    static void installPackages() {
        def packages = REQUIRED_PACKAGES.collect { "\"$it\"" }.join(",")
        execute("-e", "\"install.packages($packages)\"")
    }

    final static REQUIRED_PACKAGES = ["ggplot2", "grid", "gridExtra", "reshape", "ape", "MASS", "plotrix"]
}
