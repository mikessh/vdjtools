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
 *
 * Last modified on 7.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

public final class AaNotNtKey extends ClonotypeKey {
    public AaNotNtKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3aa().equals(other.getCdr3aa()) &&
                !clonotype.getCdr3nt().equals(other.getCdr3nt());
    }

    @Override
    public int hashCode() {
        return clonotype.getCdr3aa().hashCode();
    }
}
