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
    private class MutationsSet {
        final Collection<Mutation> mutations

        MutationsSet(Collection<Mutation> mutations) {
            this.mutations = mutations
        }

        boolean redundant = false

        int size() {
            mutations.size()
        }

        @Override
        String toString() {
            mutations.join("|")
        }
    }

    private final Map<String, Map<String, MutationsSet>> innerMap = new HashMap<>()
    private final Map<String, Integer> degreeMap = new HashMap<>()
    private final List<List<String>> subGraphs = Collections.synchronizedList(new LinkedList<>())
    final Set<Mutation> filteredShms = new HashSet<>()

    void addAll(Collection<Edge> edges) {
        synchronized (innerMap) {
            edges.each {
                addEdge(it.from, it.to, it.mutations)
            }
            //subGraphs.add([edges.collect { it.from }, edges.collect { it.to }].flatten())
            subGraphs.add(edges.collect { it.from })
        }
    }

    private void addEdge(String from, String to, Collection<Mutation> mutations) {
        def existingEdges = innerMap[from]
        if (existingEdges == null)
            innerMap.put(from, existingEdges = new HashMap<String, MutationsSet>())
        existingEdges.put(to, new MutationsSet(mutations))
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

                        if (edgeInfo != null && edgeInfo.size() == level) {
                            def connectivityCheck = new ConnectivityCheck()

                            if (connectivityCheck.checkNodes(connectivityMap, from, to)) {
                                edgeInfo.redundant = true
                                degreeMap.put(from, degreeMap[from] - edgeInfo.size())
                                degreeMap.put(to, degreeMap[to] - edgeInfo.size())
                            } else {
                                // update connectivity map
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

        //println innerMap

        filteredShms.addAll(innerMap.values().collect {
            it.values().findAll { !it.redundant }.collect { it.mutations }
        }.flatten())
    }

    Collection<EdgeInfo> collectEdges() {
        innerMap.collect { Map.Entry<String, Map<String, MutationsSet>> from ->
            from.value.findAll { to -> !to.value.redundant }.collect {
                Map.Entry<String, MutationsSet> to ->
                    to.value.mutations.collect { mutation ->
                        new EdgeInfo(from.key, to.key, mutation, to.value.size())
                    }
            }
        }.flatten()
    }
}
