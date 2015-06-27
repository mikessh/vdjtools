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
    private final int index;
    private final AtomicLong counter = new AtomicLong();
    private final Map<String, PropertyCounter> propertyCounters = new HashMap<>();

    public AminoAcidProfileBin(int index, AminoAcidPropertyGroup... aminoAcidPropertyGroups) {
        this.index = index;
        for (AminoAcidPropertyGroup aminoAcidPropertyGroup : aminoAcidPropertyGroups) {
            propertyCounters.put(aminoAcidPropertyGroup.getName(), new PropertyCounter(aminoAcidPropertyGroup));
        }
    }

    public void update(byte code, int weight) {
        counter.addAndGet(weight);

        for (PropertyCounter propertyCounter : propertyCounters.values()) {
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

        public long get(String property) {
            return counters.get(property).get();
        }

        public Map<String, Long> getAll() {
            Map<String, Long> counterMap = new HashMap<>();
            for (String property : group.getProperties()) {
                counterMap.put(property, counters.get(property).get());
            }
            return counterMap;
        }

        public AminoAcidPropertyGroup getGroup() {
            return group;
        }
    }

    public int getIndex() {
        return index;
    }

    public long getTotal() {
        return counter.get();
    }

    public Map<String, Long> getCount(String groupName) {
        PropertyCounter propertyCounter = propertyCounters.get(groupName);
        return propertyCounter.getAll();
    }

    public long getCount(String groupName, String property) {
        PropertyCounter propertyCounter = propertyCounters.get(groupName);
        return propertyCounter.get(property);
    }

    public Map<String, Map<String, Long>> getSummary() {
        Map<String, Map<String, Long>> summary = new HashMap<>();

        for (Map.Entry<String, PropertyCounter> entry : propertyCounters.entrySet()) {
            summary.put(entry.getKey(), entry.getValue().getAll());
        }
        
        return summary;
    }
}
