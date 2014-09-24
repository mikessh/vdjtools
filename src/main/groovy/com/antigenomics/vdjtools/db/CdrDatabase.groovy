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

class CdrDatabase implements Iterable<String> {
    public final String[] header
    private final HashMap<String, List<String[]>> cdrAnnotation = new HashMap<>()

    public CdrDatabase(String dbName) {
        def dbReader = CommonUtil.resourceStreamReader("db/${dbName}.txt")
        header = dbReader.readLine().split("\t")
        def line
        while ((line = dbReader.readLine()) != null) {
            def splitLine = line.split("\t")
            def cdr3aa = splitLine[0]
            def annotationList = cdrAnnotation[cdr3aa]
            if (!annotationList)
                cdrAnnotation.put(cdr3aa, annotationList = new ArrayList<String[]>())
            annotationList.add(splitLine[1..-1] as String[])
        }
    }

    public List<String[]> getAnnotationEntries(String cdr3aa) {
        cdrAnnotation[cdr3aa]
    }

    @Override
    Iterator<String> iterator() {
        cdrAnnotation.keySet().iterator()
    }
}
