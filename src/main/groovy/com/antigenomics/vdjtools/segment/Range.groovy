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
 * Last modified on 16.6.2014 by mikesh
 */

package com.antigenomics.vdjtools.segment

class Range {
    final int from, to

    Range(int from, int to) {
        this.from = from
        this.to = to
    }

    boolean contains(int x) {
        x >= from && x < to
    }

    int size() {
        to - from
    }
}
