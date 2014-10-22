package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype

/**
 * Created by mikesh on 10/22/14.
 */
class CdrDatabaseMatch {
    private static final List<MatchSubstitution> dummy = new LinkedList<>()
    private final Clonotype query
    private final boolean vMatch, jMatch
    private final List<MatchSubstitution> substitutions
    private final CdrDatabaseEntry subject

    CdrDatabaseMatch(Clonotype query, CdrDatabaseEntry subject, boolean vMatch, boolean jMatch) {
        this(query, subject, vMatch, jMatch, dummy)
    }

    CdrDatabaseMatch(Clonotype query, CdrDatabaseEntry subject, boolean vMatch, boolean jMatch,
                     List<MatchSubstitution> substitutions) {
        this.query = query
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.subject = subject
        this.substitutions = substitutions
    }

    Clonotype getQuery() {
        query
    }

    boolean getvMatch() {
        vMatch
    }

    boolean getjMatch() {
        jMatch
    }

    List<MatchSubstitution> getSubstitutions() {
        Collections.unmodifiableList(substitutions)
    }

    CdrDatabaseEntry getSubject() {
        subject
    }
}
