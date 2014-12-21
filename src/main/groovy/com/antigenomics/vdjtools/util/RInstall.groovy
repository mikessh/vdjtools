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
 * Last modified on 22.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.util

import java.util.zip.ZipInputStream

println "[RInstall] Opening resources stream"
def src = RInstall.class.protectionDomain.codeSource,
    jar = src.location,
    zip = new ZipInputStream(jar.openStream())
def entry
def dependencies = new HashSet<String>()

println "[RInstall] Scanning for dependencies"

while ((entry = zip.nextEntry)) {
    if (entry.name.toUpperCase().endsWith(".R")) {
        println "[RInstall] Scanning $entry.name"
        CommonUtil.resourceStreamReader(entry.name).readLines().each { String line ->
            if (line =~ /require\(.+\)/)
                line.split("require\\(").each { String token ->
                    if (token.contains(")")) {
                        def dependency = token.split("\\)")[0]
                        println "$dependency"
                        dependencies.add(dependency)
                    }
                }
        }
    }
}

println "[RInstall] Full list of dependencies to be installed:\n${dependencies.join(" ")}"

RUtil.install(dependencies as String[])

println "[RInstall] Testing"

RUtil.test(dependencies as String[])

println "[RInstall] Finished"