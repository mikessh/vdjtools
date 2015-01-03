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
 * Last modified on 11.12.2014 by mikesh
 */

package com.antigenomics.vdjtools.db

import com.antigenomics.vdjdb.core.db.CdrEntry
import com.antigenomics.vdjdb.core.query.CdrSearchResult
import com.antigenomics.vdjtools.Clonotype
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.AminoAcidSequence

class CdrMatch {
    private final CdrEntry subject
    private final boolean vMatch, jMatch
    private final Clonotype query
    private final Alignment<AminoAcidSequence> alignment

    public static List<CdrMatch> collectMatches(Clonotype query,
                                                CdrSearchResult cdrSearchResult,
                                                boolean requireVMatch, boolean requireJMatch) {
        requireVMatch = !requireVMatch
        requireJMatch = !requireJMatch
        def result = new ArrayList<CdrMatch>()
        cdrSearchResult.cdrEntrySet.each {
            def wrapper = new CdrMatch(it, query, cdrSearchResult.alignment)
            if ((requireVMatch || wrapper.vMatch) && (requireJMatch && wrapper.jMatch))
                result.add(wrapper)
        }
        result
    }

    public CdrMatch(CdrEntry subject, Clonotype query,
                    Alignment<AminoAcidSequence> alignment) {
        this.subject = subject
        this.query = query
        this.vMatch = query.v == subject.v //todo: re-implement V stuff
        this.jMatch = query.j == subject.j
        this.alignment = alignment
    }

    public CdrEntry getSubject() {
        subject
    }

    public Clonotype getQuery() {
        query
    }

    public boolean getvMatch() {
        vMatch
    }

    public boolean getjMatch() {
        jMatch
    }

    public Alignment<AminoAcidSequence> getAlignment() {
        alignment
    }
    public static final String HEADER = "score\tquery_cdr3aa\tquery_v\tquery_j\t" +
            "subject_cdr3aa\tsubject_v\tsubject_j\t" +
            "v_match\tj_match\tmismatches"

    @Override
    public String toString() {
        [alignment.score,
         query.cdr3aa, query.v, query.j,
         subject.cdr3aa, subject.v, subject.j,
         vMatch, jMatch, alignment.absoluteMutations.toString()
        ].flatten().join("\t")
    }
}
