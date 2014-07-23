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

enum Software {
    MiTcr("mitcr", "\t", null, 2),
    IgBlast("igblast", "\t", "#", 0),
    CdrBlast("cdrblast", "\t", "#", 0),
    Simple("simple", "\t", "#", 0)

    final String name, delimiter, comment
    final int headerLineCount

    Software(String name, String delimiter, String comment, int headerLineCount) {
        this.name = name
        this.delimiter = delimiter
        this.comment = comment
        this.headerLineCount = headerLineCount
    }

    static Software byName(String name) {
        name = name.toLowerCase()
        def software = values().find { it.name == name }
        if (!software)
            throw new IllegalArgumentException("Unknown software $name")
        software
    }
}
