package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype

/**
 * Created by mikesh on 10/22/14.
 */
class MatchSubstitution {
    public final Clonotype parent
    public final int pos
    public final char to

    MatchSubstitution(Clonotype parent, int pos, char to) {
        this.parent = parent
        this.pos = pos
        this.to = to
    }

    public char getFrom() {
        parent.cdr3aa.charAt(pos)
    }

    @Override
    String toString() {
        "$pos:$from>$to"
    }
}
