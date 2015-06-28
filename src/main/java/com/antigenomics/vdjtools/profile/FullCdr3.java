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

import com.antigenomics.vdjtools.sample.Clonotype;
import com.milaboratory.core.sequence.AminoAcidSequence;
import com.milaboratory.core.sequence.NucleotideSequence;

public class FullCdr3 implements SequenceRegion {
    @Override
    public String getName() {
        return "CDR3-full";
    }

    @Override
    public AminoAcidSequence extractAminoAcid(Clonotype clonotype) {
        return clonotype.getCdr3aaBinary();
    }

    @Override
    public NucleotideSequence extractNucleotide(Clonotype clonotype) {
        return clonotype.getCdr3ntBinary();
    }
}
