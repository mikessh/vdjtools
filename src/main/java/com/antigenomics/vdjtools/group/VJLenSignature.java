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

import com.antigenomics.vdjtools.sample.Clonotype;

public class VJLenSignature extends GroupSignature {
    public VJLenSignature(Clonotype parent) {
        super(parent);
    }

    public String getV() {
        return parent.getV();
    }

    public String getJ() {
        return parent.getJ();
    }

    public int getLength() {
        return parent.getCdr3Length();
    }

    @Override
    public int hashCode() {
        return 31 * (getV().hashCode() +
                31 * getJ().hashCode()) +
                getLength();
    }

    @Override
    public boolean equals(Object o) {
        VJLenSignature key = (VJLenSignature) o;
        return getV().equals(key.getV()) &&
                getJ().equals(key.getJ()) &&
                getLength() == key.getLength();
    }
}
