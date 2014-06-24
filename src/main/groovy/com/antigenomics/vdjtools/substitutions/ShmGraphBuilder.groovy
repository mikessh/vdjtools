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

package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.*
import groovyx.gpars.GParsPool

class ShmGraphBuilder {
    final Collection<Collection<Clonotype>> spectratype// = new HashMap<String, List<Clonotype>>()

    public ShmGraphBuilder(MutationGraph cdr3graph, ClonotypeMap clonotypeMap, double mutationRatioThreshold) {
        /*clonotypeMap.clonotypes.each { clonotype ->
            def key = clonotype.v + "\t" + clonotype.cdr3nt
            def clonotypes = spectratype[key]
            if (clonotypes == null)
                spectratype.put(key, clonotypes = new ArrayList<Clonotype>())
            clonotypes.add(clonotype)
        }*/
        spectratype = cdr3graph.getFilteredSubgraphs(mutationRatioThreshold).collect {
            new HashSet<Clonotype>(it.collect {
                clonotypeMap.getByCdr3(clonotypeMap.getByKey(it).cdr3nt)
            }.flatten())
        }
    }

    MutationGraph buildGraph() {
        def graph = new MutationGraph()

        //println "[${new Date()} INFO] Building graph"

        GParsPool.withPool Util.THREADS, {
            spectratype.eachParallel { family ->
                def edges = new LinkedList<EdgeBundle>()
                family.each { Clonotype cloneA ->
                    family.each { Clonotype cloneB ->
                        if (cloneA != cloneB) {
                            def shmIntersection = cloneA.shms.size() > cloneB.shms.size() ?
                                    cloneA.shms.intersect(cloneB.shms) : cloneB.shms.intersect(cloneA.shms)

                            def shmsBA = cloneA.shms.findAll { !shmIntersection.contains(it) }.collect { Mutation shm ->
                                shm.reassignParent(cloneB)
                            }

                            def shmsAB = cloneA.shms.findAll { !shmIntersection.contains(it) }.collect { Mutation shm ->
                                shm.reassignParent(cloneA)
                            }

                            if (shmsAB.size() > 0)
                                edges.add(new EdgeBundle(cloneA.key, cloneB.key, shmsAB, cloneA.VSegmentData.size()))

                            if (shmsBA.size() > 0)
                                edges.add(new EdgeBundle(cloneB.key, cloneA.key, shmsBA, cloneB.VSegmentData.size()))
                        }
                    }
                }
                graph.addAll(edges)
            }
        }

        //println "[${new Date()} INFO] Removing redundancy"

        graph.removeRedundancy()

        graph
    }
}
