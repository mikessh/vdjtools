package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
public abstract class ClonotypeKey {
    protected final Clonotype clonotype;

    public ClonotypeKey(Clonotype clonotype) {
        this.clonotype = clonotype;
    }

    abstract boolean equals(Clonotype other);

    @Override
    public abstract int hashCode();

    @Override
    public boolean equals(Object o) {
        // no it shouldn't
        return this.equals(((ClonotypeKey) o).clonotype);
    }
}
