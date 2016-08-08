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

package com.antigenomics.vdjtools.annotate

import com.antigenomics.vdjtools.misc.CommonUtil
import com.milaboratory.core.sequence.AminoAcidSequence

class KnownAminoAcidProperties {
    private final Map<String, AaProperty> aminoAcidPropertiesByName = new HashMap<>()

    public static final KnownAminoAcidProperties INSTANCE = new KnownAminoAcidProperties()

    private KnownAminoAcidProperties() {
        // first load simple ones
        def lines = CommonUtil.resourceStream("profile/aa_property_table.txt").readLines().findAll {
            !it.startsWith("#")
        }

        def header = lines[0].split("[\t ]+")

        if (header[0] != "amino_acid") {
            throw new RuntimeException("Simple amino acid property file header should start with 'amino_acid'.")
        }

        def propertyNames = header[1..-1]

        int n = AminoAcidSequence.ALPHABET.size()

        def propertyValues = new float[propertyNames.size()][n]

        lines[1..-1].each { line ->
            def splitLine = line.split("[\t ]+")
            byte aa = AminoAcidSequence.ALPHABET.codeFromSymbol(splitLine[0].charAt(0))

            splitLine[1..-1].eachWithIndex { String it, int ind ->
                propertyValues[ind][aa] = it.toFloat()
            }
        }

        propertyNames.eachWithIndex { String name, int ind ->
            name = name.toLowerCase()
            aminoAcidPropertiesByName.put(name,
                    new SimpleAaProperty(name, propertyValues[ind]))
        }

        // load contact model
        int span = 21, mid = span / 2i

        lines = CommonUtil.resourceStream("profile/contact_model.txt").readLines().findAll {
            !it.startsWith("#")
        }

        header = lines[0]

        if (header != "aa_tcr\tpos_norm\tprob") {
            throw new RuntimeException("Contact model should have the following header: 'aa_tcr\tpos_norm\tprob'.")
        }

        propertyValues = new float[AminoAcidSequence.ALPHABET.size()][span]

        lines[1..-1].each { line ->
            def splitLine = line.split("\t")
            byte aa = AminoAcidSequence.ALPHABET.codeFromSymbol(splitLine[0].charAt(0))
            int posNorm = splitLine[1].toInteger()
            propertyValues[aa][posNorm + mid] = splitLine[2].toFloat()
        }

        def contactEstimateProperty = new Cdr3ContactEstimate(propertyValues, 0)

        aminoAcidPropertiesByName.put(contactEstimateProperty.name.toLowerCase(), contactEstimateProperty)
    }

    AaProperty getByName(String name) {
        name = name.toLowerCase()
        if (!aminoAcidPropertiesByName.containsKey(name))
            throw new IllegalArgumentException("Bad AaProperty name '$name', allowed values: ${allowedNames}")
        aminoAcidPropertiesByName[name]
    }

    List<String> getAllowedNames() {
        aminoAcidPropertiesByName.keySet().collect()
    }

    List<AaProperty> getAll() {
        aminoAcidPropertiesByName.values().collect()
    }
}
