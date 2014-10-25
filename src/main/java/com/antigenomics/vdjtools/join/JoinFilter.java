package com.antigenomics.vdjtools.join;

/**
 * Created by mikesh on 10/25/14.
 */
public interface JoinFilter {
    public boolean pass(JointClonotype jointClonotype);
}
