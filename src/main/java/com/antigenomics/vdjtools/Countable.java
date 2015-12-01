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

package com.antigenomics.vdjtools;

/**
 * Something that has abundance data associated with it.
 */
public interface Countable {
    /**
     * Gets the number of variants associated with a given object.
     * It can be either number of convergent sub-variants associated
     * with a given clonotype or number of clonotypes in sample, or
     * the number of composite clonotypes that group several
     * clonotype sub-variants for joint/pooled samples.
     *
     * @return number of variants.
     */
    public int getDiversity();

    /**
     * Gets the number of reads associated with a given object.
     *
     * @return number of reads.
     */
    public long getCount();

    /**
     * Gets the share of reads associated with a given object.
     * Should return 1 for clonotype containers.
     *
     * @return share of reads.
     */
    public double getFreq();

    /**
     * Gets the non-normalized frequency (share of reads) for a given object.
     * For the {@link com.antigenomics.vdjtools.sample.Sample} the sum of
     * frequencies as specified in plain-text input file before normalization is returned.
     *
     * @return non-normalized share of reads.
     */
    public double getFreqAsInInput();
}
