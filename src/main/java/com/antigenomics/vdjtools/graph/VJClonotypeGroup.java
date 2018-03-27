package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.misc.Segment;

public class VJClonotypeGroup implements ClonotypeGroup {
    private final Segment v, j;

    public VJClonotypeGroup(Segment v, Segment j) {
        this.v = v;
        this.j = j;
    }

    public Segment getV() {
        return v;
    }

    public Segment getJ() {
        return j;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        VJClonotypeGroup that = (VJClonotypeGroup) o;

        if (!v.equals(that.v)) return false;
        return j.equals(that.j);
    }

    @Override
    public int hashCode() {
        int result = v.hashCode();
        result = 31 * result + j.hashCode();
        return result;
    }
}
