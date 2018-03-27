package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.sample.Clonotype;

public class VJLClonotypeGroupingFactory implements ClonotypeGroupingFactory {
    @Override
    public ClonotypeGroup getGroup(Clonotype clonotype) {
        return new VJLClonotypeGroup(clonotype.getVBinary(), clonotype.getJBinary(), clonotype.getCdr3Length());
    }
}
