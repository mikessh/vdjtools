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
    public final String[] annotationHeader
    public final String ANNOTATION_HEADER

    public CdrDatabase(String filter) {
        this(CommonUtil.resourceStreamReader("vdjdb/data/vdjdb.txt"), filter)
    }

    public CdrDatabase(String fileName, String filter) {
        this(new InputStreamReader(new FileInputStream(fileName)), filter)
    }

    public CdrDatabase(InputStreamReader dbReader, String filter) {
        def headerLine = dbReader.readLine()
        if (!headerLine.startsWith("#"))
            throw new Exception("Header line SHOULD start with '#'")
        headerLine = headerLine[1..-1]
        def header = headerLine.split("\t")

        if (filter) {
            filter.split(FILTER_MARK).each { token ->
                def columnIndex = header.findIndexOf { it.toUpperCase() == token.toUpperCase() }
                if (columnIndex >= 0)
                    filter = filter.replaceAll("$FILTER_MARK$token$FILTER_MARK", "x[$columnIndex]")
            }
        }

        def cdr3aaInd = -1, vInd = -1, jInd = -1
        def annotationHeader = new ArrayList<String>()
        def annotationIndices = new ArrayList<Integer>()

        header.eachWithIndex { String it, int ind ->
            switch (it.toUpperCase()) {
                case "CDR3AA":
                    cdr3aaInd = ind
                    break
                case "V":
                    vInd = ind
                    break
                case "J":
                    jInd = ind
                    break
                default:
                    annotationHeader.add(it)
                    annotationIndices.add(ind)
                    break
            }
        }

        if (cdr3aaInd < 0 || vInd < 0 || jInd < 0)
            throw new Exception("The following columns are MANDATORY: cdr3aa, v and j columns")

        this.annotationHeader = annotationHeader as String[]
        this.ANNOTATION_HEADER = annotationHeader.join("\t")

        def line
        while ((line = dbReader.readLine()) != null) {
            def splitLine = line.split("\t") as String[]

            if (filter && !Eval.x(splitLine, filter))
                continue

            String cdr3aa = splitLine[cdr3aaInd], v, j
            (v, j) = CommonUtil.extractVDJ(splitLine[[vInd, jInd]])

            def entryList = entriesByCdr[cdr3aa]
            if (!entryList)
                entriesByCdr.put(cdr3aa, entryList = new ArrayList<CdrDatabaseEntry>())
            entryList.add(new CdrDatabaseEntry(cdr3aa, v, j, splitLine[annotationIndices] as String[], this))
        }
    }

    public List<CdrDatabaseEntry> getAt(String cdr3aa) {
        entriesByCdr[cdr3aa] ?: []
    }
}