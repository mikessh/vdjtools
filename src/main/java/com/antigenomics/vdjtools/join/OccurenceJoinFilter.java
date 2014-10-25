package com.antigenomics.vdjtools.join;

/**
 * Created by mikesh on 10/25/14.
 */
public class OccurenceJoinFilter implements JoinFilter {
    private final int occurenceThreshold;

    public OccurenceJoinFilter() {
        this(2);
    }

    public OccurenceJoinFilter(int occurenceThreshold) {
        this.occurenceThreshold = occurenceThreshold;
    }

    public int getOccurenceThreshold() {
        return occurenceThreshold;
    }

    @Override
    public boolean pass(JointClonotype jointClonotype) {
        int detectionCounter = 0;
        for (int i = 0; i < jointClonotype.getParent().getNumberOfSamples(); i++) {
            if (jointClonotype.present(i) && ++detectionCounter == occurenceThreshold) return true;
        }
        return false;
    }
}
