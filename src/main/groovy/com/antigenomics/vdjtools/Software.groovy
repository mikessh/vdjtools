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
                                        "cdr3nt", _, _, "cdr3aa",
                                        _, "v", _, "j", _, "d",
                                        "VEnd", "DStart", "DEnd", "JStart"]),
    MiXcr("\t", null, 1, false, false, ["count", "freq",
                                        "cdr3nt", "cdr3aa",
                                        "v", "d", "j",
                                        "VEnd", "DStart", "DEnd", "JStart"]),
    MiGec("\t", null, 1, false, false, ["count", "freq",
                                        "cdr3nt", "cdr3aa",
                                        "v", "j", "d",
                                        "VEnd", "DStart", "DEnd", "JStart"]),
    HigBlast("\t", null, 1, true, false, [_, _, "count", "freq",
                                         "cdr1nt", "cdr2nt", "cdr3nt", "cdr1aa", "cdr2aa", "cdr3aa",
                                         "inFrame", "noStop", "complete",
                                         4.times { _ }].flatten()),
    ImmunoSeq("\t", null, 1, true, false, ["cdr3nt", "cdr3aa", "count", "freq", "cdr3Length",
                                           _, "v", 5.times { _ },
                                           _, "d", 5.times { _ },
                                           _, "j", 5.times { _ },
                                           6.times { _ },
                                           _, "VEnd", "DStart", "DEnd", "JStart"].flatten()),
    ImgtHighVQuest("\t", null, 1, true, true, [3.times { _ }, "v", "j", "d",
                                               9.times { _ }, "cdr3nt", 47.times { _ }, // comprehensive output
                                               "VEnd", 12.times { _ }, "DStart", "DEnd", 28.times { _ }, "JStart"
    ].flatten()),
    // Tricky for ImSEQ
    // Output is not tab-delimited table, no freq, incomplete V/J names, no D, etc
    // 'per read output' is specified to re-calculate clonotype frequencies
    ImSeq("[\t:]", null, 0, true, true, ["v", "cdr3nt", "j", "count"]),
    VDJtools("\t", "#", 0, false, false, ["count", "freq",
                                          "cdr3nt", "cdr3aa",
                                          "v", "d", "j",
                                          "VEnd", "DStart", "DEnd", "JStart"])

    static final String _ = "blank"

    final String delimiter, comment
    final boolean collapseRequired, perReadOutput
    final int headerLineCount
    @Deprecated
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
