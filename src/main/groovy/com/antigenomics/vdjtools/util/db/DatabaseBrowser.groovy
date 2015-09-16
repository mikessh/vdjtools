/*
 * Copyright (c) 2015, Bolotin Dmitry, Chudakov Dmitry, Shugay Mikhail
 * (here and after addressed as Inventors)
 * All Rights Reserved
 *
 * Permission to use, copy, modify and distribute any part of this program for
 * educational, research and non-profit purposes, by non-profit institutions
 * only, without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the following
 * three paragraphs appear in all copies.
 *
 * Those desiring to incorporate this work into commercial products or use for
 * commercial purposes should contact the Inventors using one of the following
 * email addresses: chudakovdm@mail.ru, chudakovdm@gmail.com
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE INVENTORS HAS
 * NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
 * PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY
 * PATENT, TRADEMARK OR OTHER RIGHTS.
 */

package com.antigenomics.vdjtools.util.db

import com.antigenomics.vdjdb.core.db.CdrDatabase
import com.antigenomics.vdjdb.core.query.CdrDatabaseSearcher
import com.antigenomics.vdjdb.core.query.CdrSearchResult
import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample
import com.antigenomics.vdjtools.util.ExecUtil
import com.google.common.util.concurrent.AtomicDouble
import com.milaboratory.core.tree.TreeSearchParameters
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger

/**
 *  This is a suffix tree wrapper that queries pre-defined CDR3 database 
 */
public class DatabaseBrowser {
    private final boolean vMatch, jMatch
    private final TreeSearchParameters treeSearchParameters

    /**
     * Creates a new instance of database browser with specified search parameters
     * @param vMatch if set to {@code true} Variable segments are required to match
     * @param jMatch if set to {@code true} Joining segments are required to match
     * @param treeSearchParameters suffix tree search parameters
     */
    public DatabaseBrowser(boolean vMatch, boolean jMatch,
                           TreeSearchParameters treeSearchParameters) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.treeSearchParameters = treeSearchParameters
    }

    /**
     * Creates a new instance of database browser with specified search parameters
     * @param vMatch if set to {@code true} Variable segments are required to match
     * @param jMatch if set to {@code true} Joining segments are required to match
     * @param fuzzy if set to {@code true} will allow fuzzy CDR3 match:
     *              0..2 substitutions, 0..1 insertions, 0..1 deletions but less that 3 errors in total.
     *              Will use strict CDR3 match otherwise
     */
    public DatabaseBrowser(boolean vMatch, boolean jMatch, boolean fuzzy) {
        this.vMatch = vMatch
        this.jMatch = jMatch
        this.treeSearchParameters = fuzzy ? new TreeSearchParameters(2, 1, 1, 2) : null
    }

    /**
     * Query a given sample versus a given database. Parallel implementation is used.
     * @param sample sample, all clonotypes from which will be queried
     * @param cdrDatabase database so search
     * @return a {@code BrowserResult} object
     */
    public BrowserResult query(final Sample sample, final CdrDatabase cdrDatabase) {
        def dbSearch = new CdrDatabaseSearcher(cdrDatabase, treeSearchParameters)
        def matchList = new LinkedList<CdrMatch>()
        def matchContainer = Collections.synchronizedCollection(matchList)

        // number of query clones matched and their frequency sum
        def matchDivQ = new AtomicInteger(), matchFreqQ = new AtomicDouble(),
            matchFreqQgm = new AtomicDouble()

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
                        def result = dbSearch.exact(clonotype.cdr3aa)
                        if (result)
                            found = add(clonotype, result)
                    }


                    if (found) {
                        // number of query clones matched
                        matchDivQ.incrementAndGet()

                        // match frequency
                        matchFreqQ.addAndGet(clonotype.freq)

                        // geometric mean of match frequency
                        matchFreqQgm.addAndGet(Math.log10(clonotype.freq))
                    }
                }
            }
        }

        // number of subject clones matched
        def matchDivS = matchList.collectEntries { [(it.subject): 0] }.size()
        def mathDivQVal = matchDivQ.get()

        new BrowserResult(mathDivQVal, matchDivS, matchFreqQ.get(),
                (double) (mathDivQVal > 0 ? Math.pow(10, matchFreqQgm.get() / mathDivQVal) : 0),
                matchList)
    }
}
