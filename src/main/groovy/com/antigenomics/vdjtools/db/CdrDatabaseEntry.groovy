/**
 * Created by mikesh on 10/22/14.
 */

package com.antigenomics.vdjtools.db

class CdrDatabaseEntry {
    private final String cdr3aa, v, j
    private final String[] annotation
    private final CdrDatabase parent

    CdrDatabaseEntry(String cdr3aa, String v, String j, String[] annotation, CdrDatabase parent) {
        this.cdr3aa = cdr3aa
        this.v = v
        this.j = j
        this.annotation = annotation
        this.parent = parent
    }

    String getCdr3aa() {
        return cdr3aa
    }

    String getV() {
        return v
    }

    String getJ() {
        return j
    }

    String[] getAnnotation() {
        return annotation
    }

    CdrDatabase getParent() {
        return parent
    }
}
