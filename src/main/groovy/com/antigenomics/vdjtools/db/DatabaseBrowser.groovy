package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.CommonUtil

/**
 * Created by mikesh on 10/22/14.
 */
class DatabaseBrowser {
    private final boolean vMatch, jMatch
    private final SubstitutionGenerator substitutionGenerator

    public DatabaseBrowser(boolean vMatch, boolean jMatch,
                           SubstitutionGenerator substitutionGenerator) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.substitutionGenerator = substitutionGenerator
    }

    public DatabaseBrowser(boolean vMatch, boolean jMatch) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.substitutionGenerator = new BasicSubstitutionGenerator()
    }

    Collection<CdrDatabaseMatch> query(Sample sample, CdrDatabase cdrDatabase) {
        def matchList = Collections.synchronizedCollection(new LinkedList<CdrDatabaseMatch>())

        // todo: parallel
        sample.each { Clonotype clonotype ->
            def matches = new LinkedList<CdrDatabaseMatch>()

            // find exact match
            def exactMatchEntries = cdrDatabase[clonotype.cdr3aa]
            exactMatchEntries.each {
                matches.add(new CdrDatabaseMatch(clonotype, it,
                        it.v == clonotype.v, it.j == clonotype.j))
            }

            // fuzzy matches
            if (substitutionGenerator) {
                char[] cdr3aaSeq = clonotype.cdr3aa.toCharArray()
                for (int i = 0; i < cdr3aaSeq.length; i++) {
                    char oldChar = cdr3aaSeq[i]

                    CommonUtil.AA_LIST.each { char newChar ->
                        def substitution = substitutionGenerator.create(clonotype, i, oldChar, newChar)
                        if (substitution) {
                            cdr3aaSeq[i] = newChar
                            def fuzzyMatchEntries = cdrDatabase[new String(cdr3aaSeq)]
                            fuzzyMatchEntries.each {
                                matches.add(new CdrDatabaseMatch(clonotype, it,
                                        it.v == clonotype.v, it.j == clonotype.j,
                                        [substitution]))
                            }
                        }
                    }

                    cdr3aaSeq[i] = oldChar
                }
            }

            matchList.addAll(matches.findAll {
                (!vMatch || it.vMatch) && (!jMatch || it.jMatch)
            })
        }

        matchList
    }
}
