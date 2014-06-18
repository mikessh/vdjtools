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

import com.antigenomics.vdjtools.substitutions.ConnectivityCheck

class MutationGraph {
    private class EdgeInfo {
        final Collection<Mutation> mutations

        EdgeInfo(Collection<Mutation> mutations) {
            this.mutations = mutations
        }

        boolean redundant = false

        int size() {
            mutations.size()
        }
    }

    private final Map<String, Map<String, EdgeInfo>> innerMap = new HashMap<>()
    private final Map<String, Integer> degreeMap = new HashMap<>()
    private final List<List<String>> subGraphs = new LinkedList<>()

    void addSubGraph(List<String> subGraph) {
        subGraphs.add(subGraph)
    }

    void addEdge(String from, String to, Collection<Mutation> mutations) {
        def existingEdges = innerMap[from]
        if (existingEdges == null)
            innerMap.put(from, existingEdges = new HashMap<String, EdgeInfo>())
        existingEdges.put(to, new EdgeInfo(mutations))
        degreeMap.put(from, (degreeMap[from] ?: 0) + mutations.size())
        degreeMap.put(to, (degreeMap[to] ?: 0) + mutations.size())
    }

    void removeRedundancy() {
        subGraphs.each { subGraph ->
            int n = subGraph.size(), maxLevel = 0

            subGraph.each { from ->
                innerMap[from].values().each {
                    maxLevel = Math.max(maxLevel, 2 * it.size())
                }
            }

            final def connectivityMap = new HashMap<String, List<String>>(),
                      newConnectivityMap = new HashMap<String, List<String>>()

            def addPair = { String from, String to ->
                def conList = newConnectivityMap[from]
                if (conList == null)
                    newConnectivityMap.put(from, conList = new LinkedList<String>())
                conList.add(to)
            }

            // Iteratively scan for differences by 1,2,.. mutations
            // Check connectivity graph at previous iteration to see if we're not creating cycles
            for (int level = 1; level < maxLevel; level++) {
                for (int i = 0; i < n; i++) {
                    def from = subGraph[i]
                    for (int j = i + 1; j < n; j++) {
                        def to = subGraph[j]
                        def edgeInfo = innerMap[from][to]

                        if (edgeInfo.size() == level) {
                            boolean connected = false

                            def connectivityCheck = new ConnectivityCheck()

                            if (connectivityCheck.checkNodes(connectivityMap, from, to)) {
                                edgeInfo.redundant = true
                                degreeMap.put(from, degreeMap[from] - edgeInfo.size())
                                degreeMap.put(to, degreeMap[to] - edgeInfo.size())
                            }

                            // update connectivity map
                            if (connected) {
                                addPair(from, to)
                                addPair(to, from)
                            }
                        }
                    }
                }

                // Merge connectivity map
                newConnectivityMap.each {
                    def conList = connectivityMap[it.key]
                    if (conList == null)
                        connectivityMap.put(it.key, conList = new LinkedList<String>())
                    conList.addAll(it.value)
                }
            }
        }
    }

    Collection<String> edgeStrings() {
        innerMap.collect { from -> from.value.collect { to -> ["$from (${to.value}) $to.key\t$to.value"].flatten() } }
    }
}
