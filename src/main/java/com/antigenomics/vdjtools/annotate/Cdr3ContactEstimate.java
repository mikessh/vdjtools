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

package com.antigenomics.vdjtools.annotate;

import com.milaboratory.core.sequence.AminoAcidSequence;

public class Cdr3ContactEstimate implements AaProperty {
    private final float[] valuesA, valuesB;

    public Cdr3ContactEstimate(float[] valuesA, float[] valuesB) {
        this.valuesA = valuesA;
        this.valuesB = valuesB;

        if (valuesA.length != AminoAcidSequence.ALPHABET.size()) {
            throw new IllegalArgumentException("Length of valuesA " +
                    "be equal to AminoAcidSequence.ALPHABET size.");
        }

        if (valuesB.length != AminoAcidSequence.ALPHABET.size()) {
            throw new IllegalArgumentException("Length of valuesB " +
                    "be equal to AminoAcidSequence.ALPHABET size.");
        }
    }

    @Override
    public String getName() {
        return "cdr3contact";
    }

    @Override
    public float compute(AminoAcidSequence sequence, int pos) {
        byte aa = sequence.codeAt(pos);

        float posDelta = 2.0f * ((float) pos - sequence.size() / 2.0f) / (float) sequence.size();
        posDelta *= posDelta;

        float A = valuesA[aa], B = valuesB[aa];

        return (float) Math.exp(A + B * posDelta);
    }
}
