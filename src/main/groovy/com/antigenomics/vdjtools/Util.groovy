package com.antigenomics.vdjtools
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
    static final int THREADS = Runtime.runtime.availableProcessors()

    static final String AA_LIST = "[FLSYCWPHQRIMTNKVADEGX\\*\\?]"

    static final char[] NTS = ['A', 'T', 'G', 'C']

    static String memoryFootprint() {
        final factor = 1024 * 1024 * 1024

        int maxMemory = Runtime.runtime.maxMemory() / factor,
            allocatedMemory = Runtime.runtime.totalMemory() / factor,
            freeMemory = Runtime.runtime.freeMemory() / factor

        "Memory usage: $allocatedMemory of ${maxMemory + freeMemory} GB"
    }

    static final int nt2code(char nt) {
        switch (nt) {
            case 'A':
                return 0
            case 'T':
                return 1
            case 'G':
                return 2
            case 'C':
                return 3
            default:
                return -1
        }
    }

    static final char code2nt(int code) {
        switch (code) {
            case 0:
                return 'A'
            case 1:
                return 'T'
            case 2:
                return 'G'
            case 3:
                return 'C'
            default:
                return 'N'
        }
    }

    static final char rcNt(char nt) {
        switch (nt) {
            case 'A':
                return 'T'
            case 'G':
                return 'C'
            case 'T':
                return 'A'
            case 'C':
                return 'G'
            default:
                return 'N'
        }
    }


    static String rc(String seq) {
        def chars = seq.reverse().toCharArray()
        for (int i = 0; i < chars.length; i++) {
            if (chars[i] == (char) 'A')
                chars[i] = (char) 'T'
            else if (chars[i] == (char) 'T')
                chars[i] = (char) 'A'
            else if (chars[i] == (char) 'G')
                chars[i] = (char) 'C'
            else if (chars[i] == (char) 'C')
                chars[i] = (char) 'G'
            else if (chars[i] == (char) 'N')
                chars[i] = (char) 'N'
            else if (chars[i] == (char) 'a')
                chars[i] = (char) 't'
            else if (chars[i] == (char) 't')
                chars[i] = (char) 'a'
            else if (chars[i] == (char) 'g')
                chars[i] = (char) 'c'
            else if (chars[i] == (char) 'c')
                chars[i] = (char) 'g'
            else
                chars[i] = (char) 'N'
        }
        return new String(chars)
    }

    static char codon2aa(String codon) {
        String codonUpper = codon.toUpperCase()
        switch (codonUpper) {
            case 'TTT': return 'F'
            case 'TTC': return 'F'
            case 'TTA': return 'L'
            case 'TTG': return 'L'
            case 'TCT': return 'S'
            case 'TCC': return 'S'
            case 'TCA': return 'S'
            case 'TCG': return 'S'
            case 'TAT': return 'Y'
            case 'TAC': return 'Y'
            case 'TAA': return '*'
            case 'TAG': return '*'
            case 'TGT': return 'C'
            case 'TGC': return 'C'
            case 'TGA': return '*'
            case 'TGG': return 'W'
            case 'CTT': return 'L'
            case 'CTC': return 'L'
            case 'CTA': return 'L'
            case 'CTG': return 'L'
            case 'CCT': return 'P'
            case 'CCC': return 'P'
            case 'CCA': return 'P'
            case 'CCG': return 'P'
            case 'CAT': return 'H'
            case 'CAC': return 'H'
            case 'CAA': return 'Q'
            case 'CAG': return 'Q'
            case 'CGT': return 'R'
            case 'CGC': return 'R'
            case 'CGA': return 'R'
            case 'CGG': return 'R'
            case 'ATT': return 'I'
            case 'ATC': return 'I'
            case 'ATA': return 'I'
            case 'ATG': return 'M'
            case 'ACT': return 'T'
            case 'ACC': return 'T'
            case 'ACA': return 'T'
            case 'ACG': return 'T'
            case 'AAT': return 'N'
            case 'AAC': return 'N'
            case 'AAA': return 'K'
            case 'AAG': return 'K'
            case 'AGT': return 'S'
            case 'AGC': return 'S'
            case 'AGA': return 'R'
            case 'AGG': return 'R'
            case 'GTT': return 'V'
            case 'GTC': return 'V'
            case 'GTA': return 'V'
            case 'GTG': return 'V'
            case 'GCT': return 'A'
            case 'GCC': return 'A'
            case 'GCA': return 'A'
            case 'GCG': return 'A'
            case 'GAT': return 'D'
            case 'GAC': return 'D'
            case 'GAA': return 'E'
            case 'GAG': return 'E'
            case 'GGT': return 'G'
            case 'GGC': return 'G'
            case 'GGA': return 'G'
            case 'GGG': return 'G'
            default:
                if (codonUpper.contains("N") && codonUpper.length() == 3)
                    return "X" // undefined
                else
                    return '?' // incomplete
        }
    }

    static InputStreamReader resourceStreamReader(String resourceName) {
        new InputStreamReader(Util.class.classLoader.getResourceAsStream(resourceName))
    }

    static String getSubSequence(String sequence, int from, int to) {
        String left = "", right = ""
        if (from < 0) {
            left = "N" * (-from)
            from = 0
        }
        if (to > sequence.length()) {
            right = "N" * (to - sequence.length())
            to = sequence.length()
        }

        left + sequence.substring(from, to) + right
    }
}
