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
 * Last modified on 10.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.db

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.CommonUtil
import com.google.common.util.concurrent.AtomicDouble

import java.util.concurrent.atomic.AtomicInteger


class DatabaseBrowser {
    private final boolean vMatch, jMatch
    private final SubstitutionFactory substitutionGenerator

    public DatabaseBrowser(boolean vMatch, boolean jMatch,
                           SubstitutionFactory substitutionGenerator) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.substitutionGenerator = substitutionGenerator
    }

    public DatabaseBrowser(boolean vMatch, boolean jMatch, boolean fuzzy) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.substitutionGenerator = fuzzy ? new BasicSubstitutionFactory() : null
    }

    BrowserResult query(Sample sample, CdrDatabase cdrDatabase) {
        def matchList = new LinkedList<CdrDatabaseMatch>()
        def matchContainer = Collections.synchronizedCollection(matchList)

        // number of query clones matched and their frequency sum
        def matchDivQ = new AtomicInteger(), matchFreqQ = new AtomicDouble()

        //GParsPool.withPool CommonUtil.THREADS, { // todo: parallel, for this its necessary to reformat Sample class
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

                        CommonUtil.AAS.each { char newChar ->
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

                matches = matches.findAll {
                    (!vMatch || it.vMatch) && (!jMatch || it.jMatch)
                }

                if (matches.size() > 0) {
                    matchDivQ.incrementAndGet()
                    matchFreqQ.addAndGet(clonotype.freq)
                }

                matchContainer.addAll(matches)
            }
        //}

        // number of subject clones matched
        def matchDivS = new HashSet<CdrDatabaseEntry>(matchContainer.collect { it.subject }).size()

        new BrowserResult(matchDivQ.get(), matchDivS, matchFreqQ.get(), matchList)
    }
}
