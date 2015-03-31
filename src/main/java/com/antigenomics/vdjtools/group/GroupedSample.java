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
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;

import java.util.HashMap;
import java.util.Map;

public class GroupedSample {
    private final Map<GroupSignature, Group> groups = new HashMap<>();
    private final GroupingScheme groupingScheme;
    private int total = 0;
    private final DescriptiveStatistics statistics = new DescriptiveStatistics();

    public GroupedSample(GroupingScheme groupingScheme) {
        this.groupingScheme = groupingScheme;
    }

    public void addAll(ClonotypeContainer clonotypes) {
        for (Clonotype clonotype : clonotypes) {
            add(clonotype);
        }
        summarize();
    }

    public void add(Clonotype clonotype) {
        if (clonotype.isCoding()) {
            GroupSignature groupSignature = groupingScheme.getSignature(clonotype);
            Group group = groups.get(groupSignature);
            if (group == null) {
                groups.put(groupSignature, group = new Group(groupSignature));
            }
            group.add(clonotype);
            total++;
        }
    }

    public void summarize() {
        statistics.clear();

        for (Group group : groups.values())
            statistics.addValue(Math.log10(group.size()));
    }

    public Group getGroup(Clonotype clonotype) {
        return groups.get(groupingScheme.getSignature(clonotype));
    }

    public int getTotal() {
        return total;
    }

    public boolean isRare(Group group) {
        return TestUtils.t(Math.log10(group.size()), statistics) < 0.05;
    }

    public int getNumberOfGroups() {
        return groups.size();
    }
}
