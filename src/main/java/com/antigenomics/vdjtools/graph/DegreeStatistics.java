package com.antigenomics.vdjtools.graph;

public class DegreeStatistics {
    private final int degree;
    private final long groupCount;

    public static final DegreeStatistics UNDEF = new DegreeStatistics(-1, -1);

    public DegreeStatistics(int degree, long groupCount) {
        this.degree = degree;
        this.groupCount = groupCount;
    }

    public int getDegree() {
        return degree;
    }

    public long getGroupCount() {
        return groupCount;
    }

    @Override
    public String toString() {
        return degree + "\t" + groupCount;
    }
}
