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

package com.antigenomics.vdjtools.profile

class KnownCdr3Regions implements Iterable<Cdr3Region> {
    private final Map<String, Cdr3Region> regions

    public static final KnownCdr3Regions INSTANCE = new KnownCdr3Regions()

    private KnownCdr3Regions() {
        this.regions = [new VGermline(), new DGermline(), new JGermline(),
                        new VDJunction(), new DJJunction(),
                        new VJJunction(), new FullCdr3()].collectEntries {
            [(it.name): it]
        }
    }

    public Cdr3Region getByName(String name) {
        if (!regions.containsKey(name)) {
            throw new IllegalArgumentException("[ERROR] Unknown region $name, allowed values are: ${regions.keySet()}")
        }
        regions[name]
    }

    public List<String> getRegionNames() {
        regions.keySet() as List
    }

    @Override
    Iterator<Cdr3Region> iterator() {
        regions.values().iterator()
    }
}
