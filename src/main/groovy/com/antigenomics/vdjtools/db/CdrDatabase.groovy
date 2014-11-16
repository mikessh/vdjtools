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

package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.util.CommonUtil

class CdrDatabase {
    public final String[] header
    private final HashMap<String, List<CdrDatabaseEntry>> entriesByCdr = new HashMap<>()

    public CdrDatabase(){
        this(CommonUtil.resourceStreamReader("vdjdb/data/vdjdb.txt"))
    }

    public CdrDatabase(String fileName){
        this(new InputStreamReader(new FileInputStream(fileName)))
    }

    public CdrDatabase(InputStreamReader dbReader) {
        header = dbReader.readLine().split("\t")[3..-1]
        HEADER = header.join("\t")
        def line
        while ((line = dbReader.readLine()) != null) {
            def splitLine = line.split("\t")

            String cdr3aa = splitLine[0], v, j
            (v, j) = CommonUtil.extractVDJ(splitLine[1..2])

            def entryList = entriesByCdr[cdr3aa]
            if (!entryList)
                entriesByCdr.put(cdr3aa, entryList = new ArrayList<CdrDatabaseEntry>())
            entryList.add(new CdrDatabaseEntry(cdr3aa, v, j, splitLine[3..-1] as String[], this))
        }
    }

    public final String HEADER

    public List<CdrDatabaseEntry> getAt(String cdr3aa) {
        entriesByCdr[cdr3aa] ?: []
    }
}