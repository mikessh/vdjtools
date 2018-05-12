package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.sample.Clonotype;

public class DummyClonotypeGroupingFactory implements ClonotypeGroupingFactory {
    public static final DummyClonotypeGroupingFactory INSTANCE = new DummyClonotypeGroupingFactory();

    private DummyClonotypeGroupingFactory() {

    }

    @Override
    public ClonotypeGroup getGroup(Clonotype clonotype) {
        return DummyClonotypeGroup.INSTANCE;
    }
}
