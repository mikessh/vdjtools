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

class BrowserResult implements Iterable<CdrDatabaseMatch> {
    public final int sampleDiversity, databaseDiversity
    public final double sampleFreq
    private final List<CdrDatabaseMatch> matches

    BrowserResult(int sampleDiversity, int databaseDiversity, double sampleFreq, List<CdrDatabaseMatch> matches) {
        this.sampleDiversity = sampleDiversity
        this.databaseDiversity = databaseDiversity
        this.sampleFreq = sampleFreq
        this.matches = matches
    }

    int size() {
        matches.size()
    }

    @Override
    Iterator<CdrDatabaseMatch> iterator() {
        matches.iterator()
    }

    public static
    final String HEADER = "match_size\tsample_diversity_in_matches\tdb_diversity_in_matches\tsample_freq_in_matches"

    @Override
    String toString() {
        [size(), sampleDiversity, databaseDiversity, sampleFreq].join("\t")
    }
}
