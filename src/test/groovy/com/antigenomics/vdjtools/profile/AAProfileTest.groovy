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

package com.antigenomics.vdjtools.profile

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.misc.CommonUtil
import com.milaboratory.core.sequence.AminoAcidSequence
import org.junit.Test

class AAProfileTest {
    def properties = BasicAminoAcidProperties.INSTANCE.properties

    @Test
    void test1() {
        def seq1 = CommonUtil.AAS.collect().join("")
        def profileBuilder = new AminoAcidProfile(seq1.length(), properties)

        profileBuilder.getBins().each { bin ->
            // check all properties added
            properties.each {
                assert bin.getValue(it.name) == 0
            }
        }

        profileBuilder.update(new AminoAcidSequence(seq1))

        def someProp = properties[0]

        profileBuilder.getBins().each { bin ->
            assert bin.total == 1
            assert bin.getValue(someProp.name) == someProp.getAt(seq1.charAt(bin.index))
        }
    }

    @Test
    void test2() {
        def clonotypes = [
                new Clonotype(null, 1, 1.0d,
                        [2, -1, -1, -10] as int[], "TRAV5", ".", "TRAJ48",
                        "TGTCATGAGAAATTAACCTTT",
                        "CHEKLTF", true, true, true),
                new Clonotype(null, 1, 1.0d,
                        [14, -1, -1, 19] as int[], "TRAV5", ".", "TRAJ48",
                        "TGCCTCGTGGGTGACTCGTACACGGGCAGGAGAGCACTTACTTTT",
                        "CLVGDSYTGRRALTF", true, true, true)]

        clonotypes.each { clonotype ->
            KnownCdr3Regions.INSTANCE.each { region ->
                println region.name + "\t" + region.extractAminoAcid(clonotype, true)
                println region.name + "\t" + region.extractNucleotide(clonotype)
            }
        }
    }
}
