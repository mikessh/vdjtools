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

package com.antigenomics.vdjtools.group;

import com.antigenomics.vdjtools.ClonotypeContainer;
import com.antigenomics.vdjtools.sample.Clonotype;

import java.util.HashMap;
import java.util.Map;

public class GroupedSample {
    private final Map<GroupSignature, Group> bins = new HashMap<>();
    private final GroupingScheme groupingScheme;
    private int total = 0;

    public GroupedSample(GroupingScheme groupingScheme) {
        this.groupingScheme = groupingScheme;
    }

    public void addAll(ClonotypeContainer clonotypes) {
        for (Clonotype clonotype : clonotypes) {
            add(clonotype);
        }
    }

    public void add(Clonotype clonotype) {
        GroupSignature groupSignature = groupingScheme.getSignature(clonotype);
        Group group = bins.get(groupSignature);
        if (group == null) {
            bins.put(groupSignature, group = new Group(groupSignature));
        }
        group.add(clonotype);
        total++;
    }

    public Group getBin(Clonotype clonotype) {
        return bins.get(groupingScheme.getSignature(clonotype));
    }

    public int getTotal() {
        return total;
    }

    public double getAssemblyPValue(Group group) {
        double pValue = 0;

        for (Group other : bins.values()) {
            if (other.size() > group.size()) {
                pValue++;
            } else if (other.size() == group.size()) {
                pValue += 0.5;
            }
        }

        return pValue / bins.size();
    }

    public int getNumberOfBins() {
        return bins.size();
    }
}
