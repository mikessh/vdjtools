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
