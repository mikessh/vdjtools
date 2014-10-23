/**
 * Created by mikesh on 10/22/14.
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
