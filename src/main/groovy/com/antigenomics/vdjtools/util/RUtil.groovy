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

/**
 * Class containing static utils for executing R scripts 
 */
public class RUtil {
    public static final String PACKAGES_PATH = "$ExecUtil.MY_PATH/Rpackages/" // Local R library path
    public static final String NA = "NA" // NaN in R

    /**
     * Converts a given object to numeric variable
     * @param smth object to convert
     * @return a numeric string or NA if object couldn't be converted
     */
    public static String asNumeric(smth) {
        def smthStr = smth.toString()
        if (!smthStr.isDouble())
            return NA
        def value = smthStr.toDouble()
        (value.isNaN() || value.isInfinite()) ? NA : value.toString()
    }

    /**
     * Converts a given object to logical variable
     * @param smth object to convert
     * @return "T" for true or "F" for false
     */
    public static String logical(smth) {
        smth ? "T" : "F"
    }

    /**
     * Execute a given R script 
     * @param scriptName R script name
     * @param params script parameters
     */
    public static void execute(String scriptName, String... params) {
        // Create a temp file to store the script
        def scriptRes = CommonUtil.resourceStreamReader("rscripts/$scriptName")
        scriptName = UUID.randomUUID().toString() + "_" + scriptName

        def scriptFile = new File(scriptName)

        scriptFile.withPrintWriter { pw ->
            // Set up library path correctly
            // Don't do anything if packages are not installed
            // as this would misguide R not to use /usr/ library
            if (new File(PACKAGES_PATH).exists())
                pw.println(".libPaths(\"$PACKAGES_PATH\")")

            // Write the rest of script to temp file
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

    /**
     * Install specified R packages 
     * @param dependencies names of packages to install (case-sensitive)
     */
    public static void install(String... dependencies) {
        new File(PACKAGES_PATH).mkdirs()
        execute("install.r", [PACKAGES_PATH, dependencies].flatten() as String[])
    }

    /**
     * Test if specified R packages are installed
     * @param dependencies names of packages that will be checked (case-sensitive)
     */
    public static void test(String... dependencies) {
        execute("test.r", dependencies)
    }
}
