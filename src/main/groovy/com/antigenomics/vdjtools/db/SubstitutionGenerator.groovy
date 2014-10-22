package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype

/**
 * Created by mikesh on 10/22/14.
 */
public interface SubstitutionGenerator {
    MatchSubstitution create(Clonotype parent, int pos, char from, char to)
}