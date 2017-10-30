package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.sample.Clonotype
import com.antigenomics.vdjtools.sample.Sample

class FancySpectratype {
    final double[][] spectraMatrix
    final List<Clonotype> topClonotypes
    final Spectratype spectratype

    FancySpectratype(Sample sample, int top) {
        spectratype = new Spectratype(false, false)

        topClonotypes = spectratype.addAllFancy(sample, top)

        def spectratypeHist = spectratype.getHistogram(false)

        // Prepare output table

        spectraMatrix = new double[spectratype.span][top + 1]

        for (int i = 0; i < spectratype.span; i++) {
            spectraMatrix[i][0] = spectratypeHist[i]
        }

        topClonotypes.eachWithIndex { it, ind ->
            def bin = spectratype.bin(it)
            spectraMatrix[bin][top - ind] = it.freq
        }
    }
}
