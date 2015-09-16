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

package com.antigenomics.vdjtools.diversity

/**
 * An estimate of repertoire diversity computed based on the abundances of sampled clonotypes.
 */
class DiversityEstimate {
    protected final mean, std
    protected final long numberOfReads

    public static DiversityEstimate DUMMY = new DiversityEstimate("NA", "NA", -1)

    /**
     * Creates a structure holding diversity estimate summary.
     * @param mean expected value of a diversity estimate.
     * @param std standard deviation of a diversity estimate.
     * @param numberOfReads number of reads in the sample that was analyzed.
     */
    DiversityEstimate(mean, std, long numberOfReads) {
        this.mean = mean
        this.std = std
        this.numberOfReads = numberOfReads
    }

    /**
     * Gets the mean value of diversity estimate. 
     * @return mean value.
     */
    def getMean() {
        mean
    }

    /**
     * Gets the standard deviation of diversity estimate. 
     * @return standard deviation.
     */
    def getStd() {
        std
    }

    /**
     * Gets the number of reads in the sample that was analyzed.  
     * @return associated read count.
     */
    long getNumberOfReads() {
        numberOfReads
    }

    /**
     * Plain text row for tabular output
     */
    @Override
    public String toString() {
        [mean, std].join("\t")
    }
}
