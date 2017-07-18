package com.antigenomics.vdjtools.model;

import com.antigenomics.vdjtools.sample.Clonotype;

public interface CountMatrix {
    void update(Clonotype clonotype);

    long getCount(Clonotype clonotype);
}
