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
 * Last modified on 21.6.2014 by mikesh
 */

package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.LinkedListExt

/**
 * A simple algorithm to check if two nodes are connected in graph
 */
class ConnectivityCheck {
    private final Map<String, Node> nodes = new HashMap<>()

    ConnectivityCheck(Collection<String> nodes) {
        nodes.each {
            this.nodes.put(it, new Node())
        }
    }

    void connect(String node1, String node2) {
        nodes[node1].connect(nodes[node2])
    }

    boolean connected(String node1, String node2) {
        nodes[node1].parent == nodes[node2].parent
    }

    private class Node {
        ConnectedComponent parent

        Node() {
            this.parent = new ConnectedComponent(this)
        }

        void connect(Node other) {
            if (this.parent.size() > other.parent.size())
                this.parent.merge(other.parent)
            else
                other.parent.merge(this.parent)
        }
    }

    private class ConnectedComponent {
        final LinkedListExt<Node> nodes = new LinkedListExt<>()

        ConnectedComponent(Node node) {
            nodes.add(node)
        }

        void merge(ConnectedComponent other) {
            this.nodes.concatenate(other.nodes)
            other.nodes.each { it.parent = this }
        }

        int size() {
            nodes.size()
        }
    }
}
