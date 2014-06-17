package com.antigenomics.vdjtools.substitutions

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

/**
 * A simple algorithm to check if two nodes are connected in graph
 */
class ConnectivityCheck {
    def searchedNodes = new HashSet<String>()

    // Recursively scan graph
    boolean checkNodes(HashMap<String, List<String>> graph,
                       String start, String end) {
        searchedNodes.add(start)
        def subNetwork = graph[start]

        for (String subNode : subNetwork) {
            if (!searchedNodes.contains(subNode)) {
                if (subNode == end)
                    return true
                else if (checkNodes(graph, subNode, end))
                    return true
            }
        }

        return false
    }
}
