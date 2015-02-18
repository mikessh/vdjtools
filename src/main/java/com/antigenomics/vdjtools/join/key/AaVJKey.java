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

public final class AaVJKey extends ClonotypeKey {
    public AaVJKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3aa().equals(other.getCdr3aa()) &&
                clonotype.getV().equals(other.getV()) &&
                clonotype.getJ().equals(other.getJ());
    }

    @Override
    public int hashCode() {
        return 31 * (clonotype.getCdr3aa().hashCode() * 31 + clonotype.getV().hashCode()) +
                clonotype.getJ().hashCode();
    }
}
