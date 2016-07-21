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

package com.antigenomics.vdjtools.profile;

import com.milaboratory.core.sequence.AminoAcidSequence;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import static com.milaboratory.core.sequence.AminoAcidSequence.ALPHABET;

public class AminoAcidProperty {
    public static AminoAcidProperty[] fromInput(InputStream inputStream) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));

        String line;

        do {
            line = reader.readLine();
        }
        while (line.startsWith("#"));

        String[] splitLine = line.split("\t");

        AminoAcidProperty[] groups = new AminoAcidProperty[splitLine.length - 1];

        for (int i = 0; i < groups.length; i++) {
            groups[i] = new AminoAcidProperty(splitLine[i + 1].toLowerCase());
        }

        while ((line = reader.readLine()) != null) {
            splitLine = line.split("\t");

            for (int i = 0; i < groups.length; i++) {
                groups[i].put(splitLine[0].charAt(0), Double.parseDouble(splitLine[i + 1]));
            }
        }

        return groups;
    }

    private final String name;
    private final double[] propertyByAA = new double[ALPHABET.size()];

    private AminoAcidProperty(String name) {
        this.name = name;
    }

    public void put(char aminoAcid, double value) {
        propertyByAA[ALPHABET.codeFromSymbol(aminoAcid)] = value;
    }

    public double getAt(char aminoAcid) {
        return getAt(ALPHABET.codeFromSymbol(aminoAcid));
    }

    public double getAt(byte code) {
        if (code < 0 || code >= propertyByAA.length) {
            throw new IndexOutOfBoundsException("Unknown amino acid code " + code);
        }

        return propertyByAA[code];
    }

    public double computeSum(AminoAcidSequence seq) {
        double value = 0;
        for (int i = 0; i < seq.size(); i++) {
            value += getAt(seq.codeAt(i));
        }
        return value;
    }


    public double computeSum(AminoAcidSequence seq,
                             PositionalWeighting positionalWeighting) {
        double value = 0;
        for (int i = 0; i < seq.size(); i++) {
            value += getAt(seq.codeAt(i)) *
                    positionalWeighting.getPositionWeight(i, seq.size());
        }
        return value;
    }

    public String getName() {
        return name;
    }
}
