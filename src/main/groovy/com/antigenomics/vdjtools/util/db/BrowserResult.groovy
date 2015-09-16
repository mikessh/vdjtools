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
