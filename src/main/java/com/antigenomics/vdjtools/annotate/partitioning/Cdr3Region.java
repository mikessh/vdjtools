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

package com.antigenomics.vdjtools.annotate.partitioning;

import com.antigenomics.vdjtools.sample.Clonotype;
import com.milaboratory.core.Range;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.util.Bit2Array;

public abstract class Cdr3Region implements SequenceRegion {
    private static final Range EMPTY = new Range(0, 0);

    protected abstract Range getRange(Clonotype clonotype);

    private Range safeGetRange(Clonotype clonotype) {
        Range range = getRange(clonotype);
        return (range.getFrom() >= range.getTo() || range.getFrom() < 0 || range.getTo() < 0) ? EMPTY : range;
    }

    @Override
    public AminoAcidSequence extractAminoAcid(Clonotype clonotype) {
        Range range = safeGetRange(clonotype);

        // convert to amino acids,
        // full codon belongs to a given region if at least one base of it belongs to it
        int from = range.getFrom() / 3, to = range.getTo() / 3;

        if (from >= to) {
            return new AminoAcidSequence(new byte[0]);
        }

        return clonotype.getCdr3aaBinary().getRange(new Range(from, to));
    }

    @Override
    public NucleotideSequence extractNucleotide(Clonotype clonotype) {
        Range range = safeGetRange(clonotype);

        if (range == EMPTY) {
            return new NucleotideSequence("");
        }

        return clonotype.getCdr3ntBinary().getRange(range);
    }
}
