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
 * Last modified on 23.10.2014 by mikesh
 */




package com.antigenomics.vdjtools.db

class CdrDatabaseEntry {
    public final String cdr3aa, v, j
    private final String[] annotation
    private final CdrDatabase parent

    CdrDatabaseEntry(String cdr3aa, String v, String j, String[] annotation, CdrDatabase parent) {
        this.cdr3aa = cdr3aa
        this.v = v
        this.j = j
        this.annotation = annotation
        this.parent = parent
    }

    String[] getAnnotation() {
        annotation
    }

    CdrDatabase getParent() {
        parent
    }

    boolean equals(o) {
        def that = (CdrDatabaseEntry) o
        cdr3aa == that.cdr3aa && j == that.j && v == that.v
    }

    int hashCode() {
        int result
        result = cdr3aa.hashCode()
        result = 31 * result + v.hashCode()
        31 * result + j.hashCode()
    }
}
