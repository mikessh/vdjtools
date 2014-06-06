package com.antigenomics.vdjtools

import java.util.regex.Matcher

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
class Util {
    static final String AA_LIST = "[FLSYCWPHQRIMTNKVADEGX\\*\\?]"

    static List<String> groomMatch(Matcher matcher) {
        matcher.size() > 0 ? matcher[0][1..-1] : null//[]
    }

    static String getDataPath(String dataPath) {
        if (!dataPath) {
            def SCRIPT_SOURCE = new File(getClass().protectionDomain.codeSource.location.path)
            dataPath = SCRIPT_SOURCE.parent.replaceAll("%20", " ")

            if (SCRIPT_SOURCE.absolutePath.endsWith(".groovy")) // trim /src for script
                dataPath = dataPath.replaceAll(/(?:src\/){1}.+/, "")
        } else {
            def scriptParentDir = new File(dataPath)
            if (!scriptParentDir.exists()) {
                println "Bad path to data bundle"
                System.exit(-1)
            }
            dataPath = scriptParentDir.absolutePath
        }
        dataPath
    }
}
