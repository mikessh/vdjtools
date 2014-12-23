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

import java.util.zip.GZIPInputStream

class CommonUtil {
    static final String AA_LIST = "[FLSYCWPHQRIMTNKVADEGX\\*\\?]"

    static final char[] NTS = ['A', 'T', 'G', 'C'],
                        AAS = ['F', 'L', 'S', 'Y',
                               'C', 'W', 'P', 'H',
                               'Q', 'R', 'I', 'M',
                               'T', 'N', 'K', 'V',
                               'A', 'D', 'E', 'G']

    static final byte aa2code(char aa) {
        switch (aa) {
            case 'F':
                return 0
            case 'L':
                return 1
            case 'S':
                return 2
            case 'Y':
                return 3
            case 'C':
                return 4
            case 'W':
                return 5
            case 'P':
                return 6
            case 'H':
                return 7
            case 'Q':
                return 8
            case 'R':
                return 9
            case 'I':
                return 10
            case 'M':
                return 11
            case 'T':
                return 12
            case 'N':
                return 13
            case 'K':
                return 14
            case 'V':
                return 15
            case 'A':
                return 16
            case 'D':
                return 17
            case 'E':
                return 18
            case 'G':
                return 19
            default:
                return -1
        }
    }

    static final char code2aa(byte code) {
        switch (code) {
            case 0:
                return 'F'
            case 1:
                return 'L'
            case 2:
                return 'S'
            case 3:
                return 'Y'
            case 4:
                return 'C'
            case 5:
                return 'W'
            case 6:
                return 'P'
            case 7:
                return 'H'
            case 8:
                return 'Q'
            case 9:
                return 'R'
            case 10:
                return 'I'
            case 11:
                return 'M'
            case 12:
                return 'T'
            case 13:
                return 'N'
            case 14:
                return 'K'
            case 15:
                return 'V'
            case 16:
                return 'A'
            case 17:
                return 'D'
            case 18:
                return 'E'
            case 19:
                return 'G'
            default:
                return 'X'
        }
    }

    static final byte nt2code(char nt) {
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

    static final char code2nt(byte code) {
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
        new InputStreamReader(CommonUtil.class.classLoader.getResourceAsStream(resourceName))
    }

    static InputStream getFileStream(String fileName) {
        fileName.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fileName)) : new FileInputStream(fileName)
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

    /**
     * Bring V/D/J alleles to unified look
     */
    public static List<String> extractVDJ(List<String> vdj) {
        vdj.collect {
            def major = it.split(",")[0]
            major = major.split("\\*")[0] // trim allele if present
            major = major.replaceAll("\"", "").trim()  // zap characters introduced by file editing in external software (Excel,etc)
            major.length() > 0 ? major : "."
        }
    }
}
