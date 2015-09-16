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