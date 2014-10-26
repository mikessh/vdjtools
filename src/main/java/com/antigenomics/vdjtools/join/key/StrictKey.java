package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
public final class StrictKey extends ClonotypeKey {
    public StrictKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getKey().equals(other.getKey());
    }

    @Override
    public int hashCode() {
        return clonotype.getKey().hashCode();
    }
}
