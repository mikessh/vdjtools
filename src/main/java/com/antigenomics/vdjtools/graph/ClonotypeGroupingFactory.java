package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.sample.Clonotype;

public interface ClonotypeGroupingFactory {
    ClonotypeGroup getGroup(Clonotype clonotype);
}
