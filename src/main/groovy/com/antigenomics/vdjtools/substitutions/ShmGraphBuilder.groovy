package com.antigenomics.vdjtools.substitutions

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.MutationGraph

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
class ShmGraphBuilder {
    def graph = new MutationGraph()
    def spectratype = new HashMap<String, List<Clonotype>>()

    public ShmGraphBuilder(ClonotypeMap clonotypeMap) {
        clonotypeMap.clonotypes.each { clonotype ->
            def key = clonotype.v + "\t" + clonotype.cdr3nt
            def clonotypes = spectratype[key]
            if (clonotypes == null)
                spectratype.put(key, clonotypes = new ArrayList<Clonotype>())
            clonotypes.add(clonotype)
        }
    }

    void buildGraph() {
        spectratype.values().each { family ->
            family.each { cloneA ->
                family.each { cloneB ->
                    if (cloneA != cloneB) {
                        def shmIntersection = cloneA.shms.size() > cloneB.shms.size() ?
                                cloneA.shms.intersect(cloneB.shms) : cloneB.shms.intersect(cloneA.shms)

                        def shmsBA = cloneA.shms.findAll { !shmIntersection.contains(it) }.collect {
                            it.reassignParent(cloneB)
                        },
                            shmsAB = cloneA.shms.findAll { !shmIntersection.contains(it) }.collect {
                                it.reassignParent(cloneA)
                            }

                        if (shmsAB.size() > 0)
                            graph.addEdge(cloneA.key, cloneB.key, shmsAB)

                        if (shmsBA.size() > 0)
                            graph.addEdge(cloneB.key, cloneA.key, shmsBA)
                    }
                }
            }
            graph.addSubGraph(family.collect { it.key })
        }
    }
}
