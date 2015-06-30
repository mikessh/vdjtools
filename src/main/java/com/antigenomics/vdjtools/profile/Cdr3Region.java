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

package com.antigenomics.vdjtools.profile;

import com.antigenomics.vdjtools.sample.Clonotype;
import com.milaboratory.core.Range;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.sequence.NucleotideSequence;
import com.milaboratory.util.Bit2Array;

public abstract class Cdr3Region implements SequenceRegion {
    private static final Range EMPTY = new Range(0, 0);

    protected abstract Range innerGetRange(Clonotype clonotype);

    protected Range getRange(Clonotype clonotype) {
        Range range = innerGetRange(clonotype);
        return (range.isReverse() || range.getFrom() < 0 || range.getTo() < 0) ? EMPTY : range;
    }

    @Override
    public AminoAcidSequence extractAminoAcid(Clonotype clonotype) {
        return extractAminoAcid(clonotype, false);
    }

    public AminoAcidSequence extractAminoAcid(Clonotype clonotype, boolean excludeCysPhe) {
        Range range = getRange(clonotype);

        if (range == EMPTY) {
            return new AminoAcidSequence(new byte[0]);
        }

        // convert to amino acids,
        // full codon belongs to a given region if at least one base of it belongs to it
        int from = range.getFrom() / 3, to = range.getTo() / 3;

        AminoAcidSequence cdr3aa = clonotype.getCdr3aaBinary();

        if (excludeCysPhe) {
            if (from == 0)
                from++;
            if (to == cdr3aa.size())
                to--;

            if (from >= to)
                return new AminoAcidSequence(new byte[0]);
        }

        range = new Range(from, to);

        return cdr3aa.getRange(range);
    }

    @Override
    public NucleotideSequence extractNucleotide(Clonotype clonotype) {
        Range range = getRange(clonotype);

        if (range == EMPTY) {
            return new NucleotideSequence(new Bit2Array(0));
        }

        return clonotype.getCdr3ntBinary().getRange(range);
    }
}
