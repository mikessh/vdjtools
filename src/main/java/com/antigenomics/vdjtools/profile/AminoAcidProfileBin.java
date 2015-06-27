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

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

public class AminoAcidProfileBin {
    private final AtomicLong counter = new AtomicLong();
    private final PropertyCounter[] propertyCounters;

    public AminoAcidProfileBin(AminoAcidPropertyGroup... aminoAcidPropertyGroups) {
        propertyCounters = new PropertyCounter[aminoAcidPropertyGroups.length];
        for (int i = 0; i < aminoAcidPropertyGroups.length; i++) {
            propertyCounters[i] = new PropertyCounter(aminoAcidPropertyGroups[i]);
        }
    }

    public void update(byte code, int weight) {
        counter.addAndGet(weight);

        for (PropertyCounter propertyCounter : propertyCounters) {
            propertyCounter.update(code, weight);
        }
    }

    private final class PropertyCounter {
        private final AminoAcidPropertyGroup group;
        private final Map<String, AtomicLong> counters = new HashMap<>();

        public PropertyCounter(AminoAcidPropertyGroup group) {
            this.group = group;
            for (String property : group.getProperties()) {
                counters.put(property, new AtomicLong());
            }
        }

        public void update(byte code, int weight) {
            counters.get(group.getAt(code)).addAndGet(weight);
        }

        public AminoAcidPropertyGroup getGroup() {
            return group;
        }
    }

    public long getCount() {
        return counter.get();

    }

    //public Map<String, Map<String, >>
}
