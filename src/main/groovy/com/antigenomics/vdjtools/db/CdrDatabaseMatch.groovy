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

import com.antigenomics.vdjtools.Clonotype

class CdrDatabaseMatch {
    private static final List<MatchSubstitution> dummy = new LinkedList<>()
    public final Clonotype query
    public final boolean vMatch, jMatch
    private final List<MatchSubstitution> substitutions
    public final CdrDatabaseEntry subject

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


    List<MatchSubstitution> getSubstitutions() {
        Collections.unmodifiableList(substitutions)
    }

    public static final String HEADER = "query_cdr3aa\tquery_v\tquery_j\t" +
            "subject_cdr3aa\tsubject_v\tsubject_j\t" +
            "v_match\tj_match\tsubstitutions"

    @Override
    String toString() {
        [query.cdr3aa, query.v, query.j,
         subject.cdr3aa, subject.v, subject.j,
         vMatch, jMatch, substitutions.size() > 0 ? substitutions.collect { it.toString() }.join(",") : ".",
         subject.annotation.collect()].flatten().join("\t")
    }
}
