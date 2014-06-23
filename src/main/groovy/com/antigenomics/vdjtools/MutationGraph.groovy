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
import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger

class MutationGraph {
    private class SubGraph {
        final int maxLevel
        final ConnectivityCheck connectivityCheck
        final Map<Integer, List<EdgeBundle>> edgeBundlesByLevel = new HashMap<>()
        final Map<String, Integer> degreeMap = new HashMap<>()

        SubGraph(Collection<EdgeBundle> edges) {
            connectivityCheck = new ConnectivityCheck([edges.collect { it.from },
                                                       edges.collect { it.to }].flatten())

            int maxLevel = 0

            edges.each {
                def from = it.from, to = it.to, lvl = it.size()
                degreeMap.put(from, (degreeMap[from] ?: 0) + lvl)
                degreeMap.put(to, (degreeMap[to] ?: 0) + lvl)
                def edgeList = edgeBundlesByLevel[lvl]
                if (edgeList == null)
                    edgeBundlesByLevel.put(lvl, edgeList = new LinkedList<EdgeBundle>())
                edgeList.add(it)
                maxLevel = Math.max(maxLevel, lvl)
            }
            this.maxLevel = maxLevel
        }

        void markRedundant() {
            // Iteratively scan for differences by 1,2,.. mutations
            // Check connectivity graph at previous iteration to see if we're not creating cycles
            (1..maxLevel).each { lvl ->
                def edgeList = edgeBundlesByLevel[lvl]

                if (edgeList) {
                    edgeList.each { EdgeBundle edgeBundle ->
                        if (connectivityCheck.connected(edgeBundle.from, edgeBundle.to)) {
                            edgeBundle.redundant.compareAndSet(false, true)
                            degreeMap.put(edgeBundle.from, degreeMap[edgeBundle.from] - lvl)
                            degreeMap.put(edgeBundle.to, degreeMap[edgeBundle.to] - lvl)
                        } else
                            connectivityCheck.connect(edgeBundle.from, edgeBundle.to)
                    }
                }
            }
        }
    }

    private final List<SubGraph> subGraphs = Collections.synchronizedList(new LinkedList<>())
    final Set<Mutation> filteredShms = new HashSet<>()
    final Collection<Edge> edges = new LinkedList<>()

    void addAll(Collection<EdgeBundle> edges) {
        subGraphs.add(new SubGraph(edges))
    }

    void removeRedundancy() {
        // Mark redundant edge bundles
        println "[${new Date()} GRAPH] Marking redundant edges for graph with ${subGraphs.size()} subgraphs"
        def counter = new AtomicInteger()
        GParsPool.withPool Util.THREADS, {
            subGraphs.eachParallel { SubGraph subGraph ->
                subGraph.markRedundant()
                println "[${new Date()} GRAPH] ${counter.incrementAndGet()} of ${subGraphs.size()} subgraphs processed"
            }
        }

        // Extract unique SHMs and edges for graph
        println "[${new Date()} GRAPH] Extracting non-redundant edges"
        subGraphs.each { SubGraph subGraph ->
            subGraph.edgeBundlesByLevel.values().each { edgeBundles ->
                edgeBundles.each { EdgeBundle edgeBundle ->
                    if (!edgeBundle.redundant.get()) {
                        edgeBundle.mutationSet.mutations.each { Mutation mutation ->
                            filteredShms.add(mutation)
                            edges.add(new Edge(edgeBundle.from, edgeBundle.to, mutation, edgeBundle.size()))
                        }
                    }
                }
            }
        }
    }
}
