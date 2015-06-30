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

package com.antigenomics.vdjtools.profile

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.util.CommonUtil
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
