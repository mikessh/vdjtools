package com.antigenomics.vdjtools

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
class EdgeInfo {
    final String keyFrom, keyTo, key
    final Mutation mutation
    final int weight

    EdgeInfo(String keyFrom, String keyTo, Mutation mutation, int weight) {
        this.keyFrom = keyFrom
        this.keyTo = keyTo
        this.mutation = mutation
        this.weight = weight
        this.key = keyFrom + " (" + mutation.key + ") " + keyTo
    }

    final static String EDGE_HEADER = "key\tweight\t" + Mutation.HEADER,
                        NET_HEADER = "from\tshm\tto"

    String edgeString() {
        key + "\t" + weight + "\t" + mutation
    }

    String netString() {
        keyFrom + "\t" + mutation.key + "\t" + keyTo
    }
}