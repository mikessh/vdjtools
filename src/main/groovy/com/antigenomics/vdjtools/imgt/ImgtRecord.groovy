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

class ImgtRecord {
    final String species, fullId, allele, type, functionality,
                 header, sequence

    // IMGT header
    // 0            1           2              3 4        5              6      7 8 9 10 11 12
    //>NW_001114291|TRBV24-1*01|Macaca mulatta|F|V-REGION|560072..560359|288 nt|1| | |  |  |288+39=327| | |
    //>L43137      |TRBJ1-1*01 |Macaca mulatta|F|J-REGION|749..796      |48 nt |3| | |  |  |48+0=48   | | |
    //

    public ImgtRecord(FastaRecord fastaRecord) {
        this.header = fastaRecord.header
        def splitHeader = header.split("\\|")
        this.species = splitHeader[2].replaceAll(/(\w)(\w*)/) { wholeMatch, initialLetter, restOfWord ->
            initialLetter.toUpperCase() + restOfWord
        }.replaceAll(" ", "")
        this.functionality = splitHeader[3]
        this.fullId = splitHeader[1]
        this.allele = fullId.substring(fullId.length() - 2, fullId.length())
        this.type = splitHeader[4]
        this.sequence = fastaRecord.seqeunce.toUpperCase()
    }

    @Override
    String toString() {
        ">$header\n$sequence"
    }
}