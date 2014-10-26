package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
public final class AaNotNtKey extends ClonotypeKey {
    public AaNotNtKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3aa().equals(other.getCdr3aa()) &&
                !clonotype.getCdr3nt().equals(other.getCdr3nt());
    }

    @Override
    public int hashCode() {
        return clonotype.getCdr3aa().hashCode();
    }
}
