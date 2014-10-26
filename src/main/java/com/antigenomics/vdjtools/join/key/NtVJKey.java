package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
public class NtVJKey extends ClonotypeKey {
    public NtVJKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3nt().equals(other.getCdr3nt()) &&
                clonotype.getV().equals(other.getV()) &&
                clonotype.getJ().equals(other.getJ());
    }

    @Override
    public int hashCode() {
        return 31 * (clonotype.getCdr3nt().hashCode() * 31 + clonotype.getV().hashCode()) +
                clonotype.getJ().hashCode();
    }
}
