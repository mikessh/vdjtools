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

package com.antigenomics.vdjtools.sample;

public class SequenceMatchFilter extends ClonotypeFilter {
    private final String sequence;
    private final boolean aminoAcid;
    private final int maxMismatches;

    public SequenceMatchFilter(String sequence, boolean aminoAcid, int maxMismatches) {
        super(false);
        this.sequence = sequence.toUpperCase();
        this.aminoAcid = aminoAcid;
        this.maxMismatches = maxMismatches;
    }

    @Override
    protected boolean checkPass(Clonotype clonotype) {
        String query = aminoAcid ? clonotype.getCdr3aa() : clonotype.getCdr3nt();

        if (query.length() != sequence.length()) {
            return false;
        } else if (query.equals(sequence)) {
            return true;
        }

        int nMismatches = 0;
        for (int i = 0; i < sequence.length(); i++) {
            if (sequence.charAt(i) != query.charAt(i)) {
                if (++nMismatches > maxMismatches) {
                    return false;
                }
            }
        }

        return true;
    }

    public String getSequence() {
        return sequence;
    }

    public boolean isAminoAcid() {
        return aminoAcid;
    }

    public int getMaxMismatches() {
        return maxMismatches;
    }
}
