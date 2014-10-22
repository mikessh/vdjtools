package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype

/**
 * Created by mikesh on 10/22/14.
 */
public class BasicSubstitutionGenerator implements SubstitutionGenerator {
    @Override
    public MatchSubstitution create(Clonotype parent, int pos, char from, char to) {
        (to != from &&
                pos > 0 && pos < (parent.cdr3aa.length() - 1)) ? // conserved Cys & Phe
                new MatchSubstitution(parent, pos, to) : null
    }
}