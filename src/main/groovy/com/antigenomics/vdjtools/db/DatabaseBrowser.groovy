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

import com.antigenomics.vdjdb.core.db.CdrDatabase
import com.antigenomics.vdjdb.core.query.CdrDatabaseSearcher
import com.antigenomics.vdjdb.core.query.CdrSearchResult
import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import com.google.common.util.concurrent.AtomicDouble
import com.milaboratory.core.tree.TreeSearchParameters
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger

class DatabaseBrowser {
    private final boolean vMatch, jMatch
    private final TreeSearchParameters treeSearchParameters

    public DatabaseBrowser(boolean vMatch, boolean jMatch,
                           TreeSearchParameters treeSearchParameters) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.treeSearchParameters = treeSearchParameters
    }

    public DatabaseBrowser(boolean vMatch, boolean jMatch, boolean fuzzy) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.treeSearchParameters = fuzzy ? new TreeSearchParameters(2, 1, 1, 2) : null
    }

    BrowserResult query(Sample sample, CdrDatabase cdrDatabase) {
        def dbSearch = new CdrDatabaseSearcher(cdrDatabase, treeSearchParameters)
        def matchList = new LinkedList<CdrMatch>()
        def matchContainer = Collections.synchronizedCollection(matchList)

        // number of query clones matched and their frequency sum
        def matchDivQ = new AtomicInteger(), matchFreqQ = new AtomicDouble()

        def add = { Clonotype clonotype, CdrSearchResult result ->
            matchContainer.addAll(CdrMatch.collectMatches(clonotype, result, vMatch, jMatch))
        }

        GParsPool.withPool ExecUtil.THREADS, {
            sample.eachParallel { Clonotype clonotype ->
                if (clonotype.coding) { // a mandatory internal check
                    boolean found = false

                    if (treeSearchParameters) {
                        // search with mismatches
                        dbSearch.search(clonotype.cdr3aa).each {
                            found |= add(clonotype, it)
                        }
                    } else {
                        // exact match (can also contain several associated entries)
                        found = add(clonotype, dbSearch.exact(clonotype.cdr3aa))
                    }


                    if (found) {
                        // number of query clones matched
                        matchDivQ.incrementAndGet()

                        // match frequency
                        matchFreqQ.addAndGet(clonotype.freq)
                    }
                }
            }
        }

        // number of subject clones matched
        def matchDivS = matchList.collectEntries { [(it.subject): 0] }.size()

        new BrowserResult(matchDivQ.get(), matchDivS, matchFreqQ.get(), matchList)
    }
}
