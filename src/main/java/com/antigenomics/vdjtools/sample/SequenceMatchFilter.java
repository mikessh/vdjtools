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

import org.apache.commons.lang3.StringUtils;

import java.util.regex.Pattern;

public class SequenceMatchFilter extends ClonotypeFilter {
    private final static Pattern aaAccepted = Pattern.compile("[FLSYCWPHQRIMTNKVADEGX\\*_\\[\\]]+"),
            ntAccepted = Pattern.compile("[ATGCN\\[\\]]+");

    private final String patternString;
    private final Pattern pattern;
    private final boolean aminoAcid;
    private final int maxMismatches;

    public SequenceMatchFilter(String patternString, boolean aminoAcid, int maxMismatches) {
        super(false);
        patternString = patternString.toUpperCase();

        if (!(aminoAcid ? aaAccepted : ntAccepted).matcher(patternString).matches()) {
            throw new IllegalArgumentException("Bad sequence pattern: " + patternString);
        }

        if (StringUtils.countMatches(patternString, '[') !=
                StringUtils.countMatches(patternString, ']')) {
            throw new IllegalArgumentException("Bad sequence pattern: " + patternString);
        }

        patternString = patternString.replaceAll("[XN]", ".");

        this.pattern = Pattern.compile(patternString);
        this.patternString = patternString;
        this.aminoAcid = aminoAcid;
        this.maxMismatches = maxMismatches;
    }

    @Override
    protected boolean checkPass(Clonotype clonotype) {
        String query = aminoAcid ? clonotype.getCdr3aa() : clonotype.getCdr3nt();

        return pattern.matcher(query).find();
    }

    public String getPatternString() {
        return patternString;
    }

    public boolean isAminoAcid() {
        return aminoAcid;
    }

    public int getMaxMismatches() {
        return maxMismatches;
    }
}
