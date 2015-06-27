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

import com.antigenomics.vdjtools.util.CommonUtil
import com.milaboratory.core.sequence.AminoAcidSequence
import org.junit.Test

class AAProfileTest {
    def propertyGroups = BasicAminoAcidProperties.INSTANCE.groups

    @Test
    void test1() {
        def profileBuilder = new AminoAcidProfile(10, propertyGroups)

        profileBuilder.getBins().each { bin ->
            propertyGroups.each { group ->
                group.properties.each { prop ->
                    assert bin.getCount(group.name, prop) == 0
                }
            }
        }

        def seq1 = CommonUtil.AAS.collect().join("")

        profileBuilder.update(new AminoAcidSequence(seq1))

        def propertyCounters = [:]

        profileBuilder.getBins().each { bin ->
            propertyGroups.each { group ->
                def total = 0
                group.properties.each { prop ->
                    def count = bin.getCount(group.name, prop)
                    def id = group.name + ":" + prop
                    propertyCounters.put(id, (propertyCounters[id] ?: 0) + count)
                    total += count
                }
                assert total == bin.total
            }
        }

        propertyCounters.values().each {
            assert it > 0
        }
    }
}
