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

package com.antigenomics.vdjtools.system

import com.antigenomics.vdjtools.CommonUtil

public enum RScript {
    PlotTimeSeries("plot_time_series.r")

    final String scriptName
    // todo: require, args

    RScript(String scriptName) {
        this.scriptName = scriptName
    }

    void execute(String... params) {
        // Create a temp file to store the script
        def scriptRes = CommonUtil.resourceStreamReader("rscripts/$scriptName")
        def scriptName1 = UUID.randomUUID().toString() + "_" + scriptName

        def scriptFile = new File(scriptName1)

        scriptFile.withPrintWriter { pw ->
            scriptRes.readLines().each {
                pw.println(it)
            }
        }

        scriptFile.deleteOnExit()

        // Run script
        def proc = ["Rscript", scriptName1, params].flatten().execute()

        proc.in.eachLine {
            println(it)
        }

        proc.out.close()
        proc.waitFor()

        if (proc.exitValue()) {
            println "[ERROR] ${proc.getErrorStream()}"
        }
    }
}