package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.sample.Clonotype;

public class VClonotypeGroupingFactory implements ClonotypeGroupingFactory {
    @Override
    public ClonotypeGroup getGroup(Clonotype clonotype) {
        return new VClonotypeGroup(clonotype.getVBinary());
    }
}
