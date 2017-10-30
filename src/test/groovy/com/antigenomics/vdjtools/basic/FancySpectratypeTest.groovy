package com.antigenomics.vdjtools.basic

import com.antigenomics.vdjtools.TestUtil
import org.junit.Test

class FancySpectratypeTest {
    @Test
    void test0() {
        [TestUtil.DEFAULT_SAMPLE_COLLECTION, TestUtil.SINGLE_EMPTY_SAMPLE].each { samples ->
            samples.each {
                if (it.diversity > 0) {
                    def fancySpectratype = new FancySpectratype(it, 10)

                    def values = fancySpectratype.spectraMatrix.collect { it.toList() }.flatten()

                    assert values.every { it >= 0 }
                    assert Math.abs(values.sum() - 1.0) < 1e-8
                } else {
                    // Test no error for empty
                    assert new FancySpectratype(TestUtil.SINGLE_EMPTY_SAMPLE.first(), 10).
                            spectraMatrix.collect { it.toList() }.flatten().every { it == 0 }
                }
            }
        }
    }
}
