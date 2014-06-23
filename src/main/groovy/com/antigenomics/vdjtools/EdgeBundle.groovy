/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package com.antigenomics.vdjtools

import java.util.concurrent.atomic.AtomicBoolean

class EdgeBundle {
    final MutationSet mutationSet
    final String from, to
    final AtomicBoolean redundant = new AtomicBoolean(false)

    EdgeBundle(String from, String to, Collection<Mutation> mutations) {
        this.mutationSet = new MutationSet(mutations)
        this.from = from
        this.to = to
    }

    int size() {
        mutationSet.size()
    }

    @Override
    String toString() {
        from + ">" + to + "[" + mutationSet + "]"
    }
}
