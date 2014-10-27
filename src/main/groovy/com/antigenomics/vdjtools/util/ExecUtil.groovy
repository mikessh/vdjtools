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


class ExecUtil {
    static void report(Object me, String message) { // todo: use this everywhere
        report(me, message, true)
    }

    static void report(Object me, String message, boolean verbose) { // todo: use this everywhere
        if (verbose) {
            def scriptName = me.class.canonicalName.split("\\.")[-1]
            println "[${new Date()} $scriptName] $message"
        }
    }

    static Object run(Script script, String args) {
        // perform cleanup
        def argArray = args.split(" ").
                findAll { it != " " && it != "" }.
                collect { it.replaceAll("//+", "/").toString() }
        println "Executing ${script.class.canonicalName} ${argArray.join(" ")}"
        script.binding.setVariable("args", argArray)
        script.run()
    }

    static void ensureDir(String fileName) {
        new File(fileName).absoluteFile.parentFile.mkdirs()
    }

    static String createTempDir(String outputPrefix) {
        def tmpPath = new File(outputPrefix).absolutePath + "-vdjtools-" + UUID.randomUUID().toString()

        def tmpFolderFile = new File(tmpPath)
        tmpFolderFile.mkdirs()

        //tmpFolderFile.deleteOnExit()

        tmpPath
    }
}
