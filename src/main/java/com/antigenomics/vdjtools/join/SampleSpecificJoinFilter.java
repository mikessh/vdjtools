package com.antigenomics.vdjtools.join;

/**
 * Created by mikesh on 10/25/14.
 */
public class SampleSpecificJoinFilter implements JoinFilter {
    private final int sampleIndex;

    public SampleSpecificJoinFilter() {
        this(0);
    }

    public SampleSpecificJoinFilter(int sampleIndex) {
        this.sampleIndex = sampleIndex;
    }

    public int getSampleIndex() {
        return sampleIndex;
    }

    @Override
    public boolean pass(JointClonotype jointClonotype) {
        return jointClonotype.present(sampleIndex);
    }
}
