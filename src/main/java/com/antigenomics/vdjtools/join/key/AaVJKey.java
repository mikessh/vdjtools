package com.antigenomics.vdjtools.join.key;

import com.antigenomics.vdjtools.Clonotype;

/**
 * Created by mikesh on 10/26/14.
 */
public final class AaVJKey extends ClonotypeKey {
    public AaVJKey(Clonotype clonotype) {
        super(clonotype);
    }

    @Override
    public boolean equals(Clonotype other) {
        return clonotype.getCdr3aa().equals(other.getCdr3aa()) &&
                clonotype.getV().equals(other.getV()) &&
                clonotype.getJ().equals(other.getJ());
    }

    @Override
    public int hashCode() {
        return 31 * (clonotype.getCdr3aa().hashCode() * 31 + clonotype.getV().hashCode()) +
                clonotype.getJ().hashCode();
    }
}
