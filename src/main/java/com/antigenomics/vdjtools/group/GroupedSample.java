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
