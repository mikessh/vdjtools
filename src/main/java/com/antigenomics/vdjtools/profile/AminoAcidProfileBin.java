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

import com.google.common.util.concurrent.AtomicDouble;

import java.util.HashMap;
import java.util.Map;

public class AminoAcidProfileBin {
    private final int index;
    private final AtomicDouble counter = new AtomicDouble();
    private final Map<String, PropertyCounter> propertyCounters = new HashMap<>();

    public AminoAcidProfileBin(int index, AminoAcidProperty... aminoAcidProperties) {
        this.index = index;
        for (AminoAcidProperty aminoAcidProperty : aminoAcidProperties) {
            propertyCounters.put(aminoAcidProperty.getName(), new PropertyCounter(aminoAcidProperty));
        }
    }

    protected void update(byte code, double weight) {
        counter.addAndGet(weight);

        for (PropertyCounter propertyCounter : propertyCounters.values()) {
            propertyCounter.update(code, weight);
        }
    }

    protected final class PropertyCounter {
        private final AminoAcidProperty group;
        private final AtomicDouble counter = new AtomicDouble();

        public PropertyCounter(AminoAcidProperty group) {
            this.group = group;
        }

        public void update(byte code, double weight) {
            counter.addAndGet(group.getAt(code) * weight);
        }

        public double getValue() {
            return counter.get();
        }

        public AminoAcidProperty getGroup() {
            return group;
        }
    }

    public int getIndex() {
        return index;
    }

    public double getTotal() {
        return counter.get();
    }

    public double getValue(String groupName) {
        PropertyCounter propertyCounter = propertyCounters.get(groupName);
        return propertyCounter.getValue();
    }

    public Map<String, Double> getSummary() {
        Map<String, Double> summary = new HashMap<>();

        for (PropertyCounter propertyCounter : propertyCounters.values()) {
            summary.put(propertyCounter.group.getName(),
                    propertyCounter.getValue());
        }

        return summary;
    }
}
