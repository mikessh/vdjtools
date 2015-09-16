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

package com.antigenomics.vdjtools.sample

/**
 * An object holding some basic statistics for a sample collection 
 */
public class SampleStatistics {
    private final long minCount, maxCount
    private final double minFreq, maxFreq
    private final int minDiversity, maxDiversity

    /**
     * Initializes a new sample statistics object with pre-computed values 
     * @param minCount minimal number of reads in a sample from given sample collection
     * @param maxCount maximal number of reads in a sample from given sample collection
     * @param minFreq minimal total frequency of clonotypes in a sample from given sample collection
     * @param maxFreq maximal total frequency of clonotypes in a sample from given sample collection
     * @param minDiversity minimal number of clonotypes in a sample from given sample collection
     * @param maxDiversity maximal number of clonotypes in a sample from given sample collection
     */
    public SampleStatistics(long minCount, long maxCount,
                            double minFreq, double maxFreq,
                            int minDiversity, int maxDiversity) {
        this.minCount = minCount
        this.maxCount = maxCount
        this.minFreq = minFreq
        this.maxFreq = maxFreq
        this.minDiversity = minDiversity
        this.maxDiversity = maxDiversity
    }

    /**
     * Gets the minimal number of reads in a sample from given sample collection
     * @return size of the smallest sample in sample collection
     */
    public long getMinCount() {
        minCount
    }

    /**
     * Gets the maximal number of reads in a sample from given sample collection
     * @return size of the largest sample in sample collection
     */
    public long getMaxCount() {
        maxCount
    }

    public double getMinFreq() {
        minFreq
    }

    public double getMaxFreq() {
        maxFreq
    }

    public int getMinDiversity() {
        minDiversity
    }

    public int getMaxDiversity() {
        maxDiversity
    }
}
