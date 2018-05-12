package com.antigenomics.vdjtools.graph;

public class DegreeStatistics {
    private final int degree;
    private final long primaryGroupCount, secondaryGroupCount;

    public static final DegreeStatistics UNDEF = new DegreeStatistics(-1, -1, -1);

    public DegreeStatistics(int degree, long primaryGroupCount, long secondaryGroupCount) {
        this.degree = degree;
        this.primaryGroupCount = primaryGroupCount;
        this.secondaryGroupCount = secondaryGroupCount;
    }

    public int getDegree() {
        return degree;
    }

    public long getPrimaryGroupCount() {
        return primaryGroupCount;
    }

    public long getSecondaryGroupCount() {
        return secondaryGroupCount;
    }

    @Override
    public String toString() {
        return degree + "\t" + primaryGroupCount + "\t" + secondaryGroupCount;
    }
}
