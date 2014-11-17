/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Last modified on 17.11.2014 by mikesh
 */



package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.util.CommonUtil

class CdrDatabase {
    public static final String FILTER_MARK = "__"
    private final HashMap<String, List<CdrDatabaseEntry>> entriesByCdr = new HashMap<>()
    public final String[] header
    public final String HEADER

    public CdrDatabase(String filter) {
        this(CommonUtil.resourceStreamReader("vdjdb/data/vdjdb.txt"), filter)
    }

    public CdrDatabase(String fileName, String filter) {
        this(new InputStreamReader(new FileInputStream(fileName)), filter)
    }

    public CdrDatabase(InputStreamReader dbReader, String filter) {
        def header = dbReader.readLine().split("\t")

        if (filter) {
            filter.split(FILTER_MARK).each { token ->
                println token
                def columnIndex = header.findIndexOf { it.toUpperCase() == token.toUpperCase() }
                if (columnIndex >= 0)
                    filter = filter.replaceAll("$FILTER_MARK$token$FILTER_MARK", "x[$columnIndex]")
            }
        }

        this.header = header[3..-1]
        this.HEADER = header.join("\t")

        def line
        while ((line = dbReader.readLine()) != null) {
            def splitLine = line.split("\t") as String[]

            if (filter && !Eval.x(splitLine, filter))
                continue

            // todo: more flexible v/j/cdr3 column search
            String cdr3aa = splitLine[0], v, j
            (v, j) = CommonUtil.extractVDJ(splitLine[1..2])

            def entryList = entriesByCdr[cdr3aa]
            if (!entryList)
                entriesByCdr.put(cdr3aa, entryList = new ArrayList<CdrDatabaseEntry>())
            entryList.add(new CdrDatabaseEntry(cdr3aa, v, j, splitLine[3..-1] as String[], this))
        }
    }

    public List<CdrDatabaseEntry> getAt(String cdr3aa) {
        entriesByCdr[cdr3aa] ?: []
    }
}