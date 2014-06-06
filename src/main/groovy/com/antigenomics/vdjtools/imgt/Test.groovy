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

package com.antigenomics.vdjtools.imgt

import com.antigenomics.vdjtools.io.FastaRecord

def record1 = new ImgtRecord(new FastaRecord(
        "J00256|IGHJ1*01|Homo sapiens|F|J-REGION|723..774|52 nt|1| | | | |52+0=52| | |",
        "" +
                "gctgaatacttccagcactggggccagggcaccctggtcaccgtctcctcag"
)
)
def record2 = new ImgtRecord(new FastaRecord(
        "M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320| | |",
        "" +
                "caggttcagctggtgcagtctggagct...gaggtgaagaagcctggggcctcagtgaag" +
                "gtctcctgcaaggcttctggttacaccttt............accagctatggtatcagc" +
                "tgggtgcgacaggcccctggacaagggcttgagtggatgggatggatcagcgcttac..." +
                "...aatggtaacacaaactatgcacagaagctccag...ggcagagtcaccatgaccaca" +
                "gacacatccacgagcacagcctacatggagctgaggagcctgagatctgacgacacggcc" +
                "gtgtattactgtgcgagaga"
)
)

def record3 = new ImgtRecord(new FastaRecord(
        "X62111|IGHV2-5*01|Homo sapiens|F|V-REGION|214..514|301 nt|1| | | | |301+21=322| | |",
        "" +
                "cagatcaccttgaaggagtctggtcct...acgctggtgaaacccacacagaccctcacg" +
                "ctgacctgcaccttctctgggttctcactcagc......actagtggagtgggtgtgggc" +
                "tggatccgtcagcccccaggaaaggccctggagtggcttgcactcatttattggaat..." +
                "......gatgataagcgctacagcccatctctgaag...agcaggctcaccatcaccaag" +
                "gacacctccaaaaaccaggtggtccttacaatgaccaacatggaccctgtggacacagcc" +
                "acatattactgtgcacacagac"
)
)

def record4 = new ImgtRecord(new FastaRecord(
        "AE000659|TRAV10*01|Homo sapiens|F|V-REGION|81133..81412|280 nt|1| | | | |280+42=322| | |",
        "" +
                "aaaaaccaagtggagcagagtcctcagtccctgatcatcctggagggaaagaactgcact" +
                "cttcaatgcaattatacagtgagcccc..................ttcagcaacttaagg" +
                "tggtataagcaagatactgggagaggtcctgtttccctgacaatcatgactttcagt..." +
                "......gagaacacaaagtcgaac...............ggaagatatacagcaactctg" +
                "gatgcagacacaaagcaaagctctctgcacatcacagcctcccagctcagcgattcagcc" +
                "tcctacatctgtgtggtgagcg"
)
)

int IMGT_V_REF = 312
int imgtGaps = record4.sequence.substring(0, IMGT_V_REF).count(".")
int MIGEC_V_REF = IMGT_V_REF - imgtGaps
println MIGEC_V_REF
println record4.sequence.substring(IMGT_V_REF - 3, IMGT_V_REF)

println "\n\n"

println ImgtToMigecParser.getJReferencePoint(record1)
