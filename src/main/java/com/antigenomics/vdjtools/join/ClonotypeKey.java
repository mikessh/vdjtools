package com.antigenomics.vdjtools.join;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
abstract class ClonotypeKey {
    protected final Clonotype clonotype;

    protected ClonotypeKey(Clonotype clonotype) {
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
