package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
public final class NtVKey extends ClonotypeKey {
    public NtVKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3nt().equals(other.getCdr3nt()) &&
                clonotype.getV().equals(other.getV());
    }

    @Override
    public int hashCode() {
        return clonotype.getCdr3nt().hashCode() * 31 + clonotype.getV().hashCode();
    }
}