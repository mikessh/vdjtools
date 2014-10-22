package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype

/**
 * Created by mikesh on 10/22/14.
 */
class MatchSubstitution {
    private final Clonotype parent
    public final int pos
    public final char to

    MatchSubstitution(Clonotype parent, int pos, char to) {
        this.parent = parent
        this.pos = pos
        this.to = to
    }

    Clonotype getParent() {
        parent
    }
}
