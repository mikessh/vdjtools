package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.misc.Segment;

public class VClonotypeGroup implements ClonotypeGroup {
    private final Segment v;

    public VClonotypeGroup(Segment v) {
        this.v = v;
    }

    public Segment getV() {
        return v;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VClonotypeGroup that = (VClonotypeGroup) o;

        return v.equals(that.v);
    }

    @Override
    public int hashCode() {
        return v.hashCode();
    }
}
