package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
public final class AaKey extends ClonotypeKey {
    public AaKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3aa().equals(other.getCdr3aa());
    }

    @Override
    public int hashCode() {
        return clonotype.getCdr3aa().hashCode();
    }
}