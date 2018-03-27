package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.sample.Clonotype;

public class VJClonotypeGroupingFactory implements ClonotypeGroupingFactory {
    @Override
    public ClonotypeGroup getGroup(Clonotype clonotype) {
        return new VJClonotypeGroup(clonotype.getVBinary(), clonotype.getJBinary());
    }
}
