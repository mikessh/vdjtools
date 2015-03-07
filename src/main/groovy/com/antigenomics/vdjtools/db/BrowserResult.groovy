/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 */

package com.antigenomics.vdjtools.db

/**
 * A class representing the result of querying CDR3 database with a list of clonotypes.
 *
 * Contains summary statistics:
 * @{code match_size}             Total number of matches between sample and database.
 *                                If a given clonotype has N matches, it will be counted N times.
 * @{sample_diversity_in_matches} Number of unique clonotypes in the sample that have at least one match
 * @{db_diversity_in_matches}     Number of entries in the database that were matched at least once
 * @{sample_freq_in_matches}      Total frequency of unique clonotypes in the sample that have at least one match
 * @{mean_matched_clone_size}     Geometric mean of frequency of unique clonotypes in the sample that have at least one match
 *
 * Using this class one can also iterate over the list of matches
 */
public class BrowserResult implements Iterable<CdrMatch> {
    private final int sampleDiversity, databaseDiversity
    private final double sampleFreq, meanCloneSize
    private final List<CdrMatch> matches

    /**
     * Creates CDR3 database batch query result
     * @param sampleDiversity Number of unique clonotypes in the sample that have at least one match
     * @param databaseDiversity Number of entries in the database that were matched at least once
     * @param sampleFreq Total frequency of unique clonotypes in the sample that have at least one match
     * @param meanCloneSize Geometric mean of frequency of unique clonotypes in the sample that have at least one match
     * @param matches list of matches
     */
    public BrowserResult(int sampleDiversity, int databaseDiversity,
                         double sampleFreq, double meanCloneSize,
                         List<CdrMatch> matches) {
        this.sampleDiversity = sampleDiversity
        this.databaseDiversity = databaseDiversity
        this.sampleFreq = sampleFreq
        this.meanCloneSize = meanCloneSize
        this.matches = matches
    }

    /**
     * Gets the number of matches in the result 
     * @return
     */
    public int size() {
        matches.size()
    }

    @Override
    /**
     * Iterate over matches 
     */
    public Iterator<CdrMatch> iterator() {
        matches.iterator()
    }

    /**
     * Gets the number of unique clonotypes in the sample that have at least one match 
     * @return
     */
    public int getSampleDiversity() {
        return sampleDiversity
    }

    /**
     * Gets the number of entries in the database that were matched at least once 
     * @return
     */
    public int getDatabaseDiversity() {
        return databaseDiversity
    }

    /**
     * Gets the total frequency of unique clonotypes in the sample that have at least one match 
     * @return
     */
    public double getSampleFreq() {
        return sampleFreq
    }

    /**
     * Gets the geometric mean of frequency of unique clonotypes in the sample that have at least one match 
     * @return
     */
    double getMeanCloneSize() {
        return meanCloneSize
    }

    /**
     * Header string, used for tabular output
     */
    public static final String HEADER = "match_size\tsample_diversity_in_matches\tdb_diversity_in_matches\t" +
            "sample_freq_in_matches\tmean_matched_clone_size"

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        [size(), sampleDiversity, databaseDiversity, sampleFreq, meanCloneSize].join("\t")
    }
}
