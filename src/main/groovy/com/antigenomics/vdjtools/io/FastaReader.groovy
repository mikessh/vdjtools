package com.antigenomics.vdjtools.io
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
class FastaReader implements Iterable<FastaRecord> {
    final BufferedReader reader
    String header = ""

    FastaReader(String fileName) {
        this(new FileReader(fileName))
    }

    FastaReader(InputStreamReader input) {
        reader = new BufferedReader(input)
    }

    FastaRecord next() {
        if (header == "") // first sequence
            header = reader.readLine()

        if (!header) // EOF
            return null

        def seq = ""
        while (true) {
            def line = reader.readLine()
            if (line == null) {
                // EOF next read, return last sequence
                def _header = header
                header = null
                return new FastaRecord(_header, seq)
            } else if (line.startsWith(">")) {
                // reset header and return current read
                def _header = header
                header = line
                return new FastaRecord(_header, seq)
            } else // sequence line
                seq += line
        }


    }

    Iterator iterator() {
        return [ hasNext: { header != null}, next: { next() } ] as Iterator
    }
}
