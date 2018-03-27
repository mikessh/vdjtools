package com.antigenomics.vdjtools.graph;

import com.antigenomics.vdjtools.misc.Segment;

public class VJLClonotypeGroup extends VJClonotypeGroup {
    private final int length;

    public VJLClonotypeGroup(Segment v, Segment j, int length) {
        super(v, j);
        this.length = length;
    }

    public int getLength() {
        return length;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        VJLClonotypeGroup that = (VJLClonotypeGroup) o;

        return length == that.length;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + length;
        return result;
    }
}
