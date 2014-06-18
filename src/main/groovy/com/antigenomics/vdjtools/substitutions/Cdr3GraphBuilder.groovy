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

import com.antigenomics.vdjtools.Clonotype
import com.antigenomics.vdjtools.Mutation
import com.antigenomics.vdjtools.MutationGraph
import com.antigenomics.vdjtools.Util

class Cdr3GraphBuilder {
    def spectratype = new HashMap<String, Map<String, List<Clonotype>>>()
    def graph = new MutationGraph()

    public Cdr3GraphBuilder(ClonotypeMap clonotypeMap) {
        clonotypeMap.clonotypes.each { clonotype ->
            def key = clonotype.v + "\t" + clonotype.cdr3nt.length()
            def clonotypesByCdr3 = spectratype[key]
            if (clonotypesByCdr3 == null)
                spectratype.put(key, clonotypesByCdr3 = new HashMap<String, List<Clonotype>>())
            def clonotypeList = clonotypesByCdr3[clonotype.cdr3nt]
            if (clonotypeList == null)
                clonotypesByCdr3.put(clonotype.cdr3nt, clonotypeList = new LinkedList<Clonotype>())
            clonotypeList.add(clonotype)
        }
    }

    void buildGraph() {
        spectratype.values().each { spectraPeak ->
            def subGraph = new ArrayList<String>()
            def cloneLists = spectraPeak.values()
            for (int i = 0; i < cloneLists.size(); i++) {
                def clonesA = cloneLists[i]
                for (int j = i + 1; j < cloneLists.size(); j++) {
                    def clonesB = cloneLists[j]
                    // select which clonotypes to connect based on SHM intersection
                    Clonotype cloneA, cloneB
                    (cloneA, cloneB) = [clonesA, clonesB].combinations().max {
                        Clonotype cl1 = it[0], cl2 = it[1]

                        cl1.shms.size() > cl2.shms.size() ?
                                cl1.shms.intersect(cl2.shms).size() :
                                cl2.shms.intersect(cl1.shms).size()
                    }.collect()
                    def mutations = extractMutations(cloneA, cloneB)
                    if (mutations != null) {
                        subGraph.add(cloneA.key)
                        subGraph.add(cloneB.key)
                        graph.addEdge(cloneA.key, cloneB.key, extractMutations(cloneB, cloneA))
                    }
                }
            }
            graph.addSubGraph(subGraph)
        }
    }

    private static Set<Mutation> extractMutations(Clonotype clone1, Clonotype clone2) {
        def mutations = new HashSet<Mutation>()
        int len = clone1.cdr3nt.length(), depth = scanDept(len)

        def mutationPositions = new LinkedList<Integer>()
        for (int i = 0; i < len; i++) {
            if (clone1.cdr3nt.charAt(i) != clone2.cdr3nt.charAt(i)) {
                mutationPositions.add(i)
                if (mutationPositions.size() > depth)
                    return null
            }
        }

        mutationPositions.each { int pos ->
            int codonStart = pos - pos % 3
            if (len < codonStart + 3) {
                String fromCodon = clone1.cdr3nt.substring(codonStart, codonStart + 3),
                       toCodon = clone2.cdr3nt.substring(codonStart, codonStart + 3)

                char fromAA = Util.codon2aa(fromCodon), toAA = Util.codon2aa(toCodon)

                mutations.add(Mutation.cdr3Mutation(pos,
                        clone1.cdr3nt.charAt(pos), clone2.cdr3nt.charAt(pos),
                        (int) (pos / 3), fromAA, toAA, false, clone1))
            }
        }
    }

    static int scanDept(int len) {
        return 5
    }
}
