package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.sample.Clonotype;

public class DummyClonotypeGroupingFactory implements ClonotypeGroupingFactory {
    @Override
    public ClonotypeGroup getGroup(Clonotype clonotype) {
        return DummyClonotypeGroup.INSTANCE;
    }
}
