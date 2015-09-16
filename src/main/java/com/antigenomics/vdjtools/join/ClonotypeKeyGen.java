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

package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.join.key.*;
import com.antigenomics.vdjtools.overlap.OverlapType;
import com.antigenomics.vdjtools.sample.Clonotype;
import com.antigenomics.vdjtools.sample.Sample;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.HashSet;
import java.util.Set;

public class ClonotypeKeyGen {
    private final OverlapType overlapType;

    public ClonotypeKeyGen() {
        this(OverlapType.Strict);
    }

    public ClonotypeKeyGen(OverlapType overlapType) {
        this.overlapType = overlapType;
    }

    public Set<ClonotypeKey> generateKeySet(Sample sample) {
        Set<ClonotypeKey> keySet = new HashSet<>();
        for (Clonotype clonotype : sample) {
            keySet.add(generateKey(clonotype));
        }
        return keySet;
    }

    public ClonotypeKey generateKey(Clonotype clonotype) {
        switch (overlapType) {
            case Nucleotide:
                return new NtKey(clonotype);

            case NucleotideV:
                return new NtVKey(clonotype);

            case NucleotideVJ:
                return new NtVJKey(clonotype);

            case AminoAcid:
                return new AaKey(clonotype);

            case AminoAcidV:
                return new AaVKey(clonotype);

            case AminoAcidVJ:
                return new AaVJKey(clonotype);

            case AminoAcidNonNucleotide:
                return new AaNotNtKey(clonotype);

            case Strict:
                return new StrictKey(clonotype);

            default:
                throw new NotImplementedException();
        }
    }
}
