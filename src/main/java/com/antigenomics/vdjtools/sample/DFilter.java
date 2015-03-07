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

import com.antigenomics.vdjtools.Clonotype;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

public class DFilter extends ClonotypeFilter {
    private final Set<String> dSegmentSet = new HashSet<>();

    public DFilter(boolean negative, String... dSegment) {
        super(negative);
        Collections.addAll(this.dSegmentSet, dSegment);
    }

    public DFilter(String... dSegment) {
        this(false, dSegment);
    }

    @Override
    protected boolean checkPass(Clonotype clonotype) {
        return dSegmentSet.contains(clonotype.getD());
    }

    public Set<String> getdSegmentSet() {
        return Collections.unmodifiableSet(dSegmentSet);
    }
}
