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

package com.antigenomics.vdjtools

enum Software {
    // todo : output mutations
    MiTcr("\t", null, 2, false, false, ["count", "freq",
                                        "cdr3nt", "blank", "blank", "cdr3aa",
                                        "blank", "v", "blank", "j", "blank", "d",
                                        "VEnd", "DStart", "DEnd", "JStart"]),
    MiGec("\t", null, 1, false, false, ["count", "freq",
                                        "cdr3nt", "cdr3aa",
                                        "v", "j", "d",
                                        "VEnd", "DStart", "DEnd", "JStart"]),
    IgBlast("\t", "#", 0, false, false, ["blank", "blank", "count", "freq",
                                         "cdr1nt", "cdr2nt", "cdr3nt", "cdr1aa", "cdr2aa", "cdr3aa",
                                         "inFrame", "noStop", "complete",
                                         "blank", "blank", "blank",
                                         "blank"]),
    ImmunoSeq("\t", null, 1, true, false, ["cdr3nt", "cdr3aa", "count", "freq", "cdr3Length",
                                           "blank", "v", "blank", "blank", "blank", "blank", "blank",
                                           "blank", "d", "blank", "blank", "blank", "blank", "blank",
                                           "blank", "j", "blank", "blank", "blank", "blank", "blank",
                                           "blank", "blank", "blank", "blank", "blank", "blank",
                                           "blank", "VEnd", "DStart", "DEnd", "JStart"]),
    Simple("\t", "#", 0, false, false, ["count", "freq",
                                        "cdr3nt", "cdr3aa",
                                        "v", "d", "j",
                                        "VEnd", "DStart", "DEnd", "JStart"])

    final String delimiter, comment
    final boolean collapseRequired, perReadOutput
    final int headerLineCount
    final List<String> printFields

    Software(String delimiter, String comment, int headerLineCount, boolean collapseRequired, boolean perReadOutput, List<String> printFields) {
        this.delimiter = delimiter
        this.comment = comment
        this.headerLineCount = headerLineCount
        this.collapseRequired = collapseRequired
        this.perReadOutput = perReadOutput
        this.printFields = printFields

        if (perReadOutput && !collapseRequired)
            throw new RuntimeException("Collapsing should always be performed for a software with per-read output")
    }

    static Software byName(String name) {
        def software = values().find {
            it.toString().toLowerCase() == name.toLowerCase()
        }
        if (!software)
            throw new IllegalArgumentException("Unknown software $name")
        software
    }

    static String allowedNames = values().collect { it.toString().toLowerCase() }.join(",")
}
