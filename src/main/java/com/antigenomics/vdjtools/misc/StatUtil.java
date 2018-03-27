package com.antigenomics.vdjtools.misc;

import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import static org.apache.commons.math3.special.Beta.regularizedBeta;
import static org.apache.commons.math3.util.CombinatoricsUtils.binomialCoefficientLog;

/**
 * Basic distributions (PDF/CDF), copied from MAGERI
 */
public final class StatUtil {
    private StatUtil() {

    }

    public static double binomialPdf(int k, int n, double p) {
        return Math.exp(binomialCoefficientLog(n, k) +
                (n - k) * Math.log(1.0 - p) + k * Math.log(p));
    }

    public static double binomialCdf(int k, int n, double p) {
        return 1.0 - regularizedBeta(1 - p, n - k, k + 1);
    }

    public static double binomialPValue(int k, int n, double p) {
        return binomialCdf(k, n, p) + 0.5 * binomialPdf(k, n, p);
    }

    public static double binomialPValue2(int k, int n, double p) {
        return Math.min(binomialPValue(k, n, p), 1 - binomialPValue(k, n, p));
    }

    /**
     * Computes the PDF of the distribution of number of successes in a sequence of
     * iid Bernoulli trials before a specified number of failures occurs.
     *
     * @param k number of successes
     * @param r number of failures until the experiment is stopped
     * @param p success probability in each experiment
     * @return probability
     */
    public static double negativeBinomialPdf(int k, int r, double p) {
        return Math.exp(binomialCoefficientLog(k + r - 1, k) +
                r * Math.log(1.0 - p) + k * Math.log(p));
    }

    /**
     * Computes the CDF of the distribution of number of successes in a sequence of
     * iid Bernoulli trials before a specified number of failures occurs.
     *
     * @param k number of successes
     * @param r number of failures until the experiment is stopped
     * @param p success probability in each experiment
     * @return cumulative
     */
    public static double negativeBinomialCdf(int k, double r, double p) {
        return 1.0 - regularizedBeta(p, k + 1, r);
    }

    public static double betaBinomialPdf(int k, int n, double alpha, double beta) {
        return Math.exp(
                (Beta.logBeta(k + alpha, n - k + beta) - Beta.logBeta(alpha, beta)) +
                        (Gamma.logGamma(n + 1) - Gamma.logGamma(k + 1) - Gamma.logGamma(n - k + 1))
        );
    }

    public static double betaBinomialCdf(int k, int n, double alpha, double beta) {
        double sum = 0;

        for (int i = 0; i < k; i++) {
            sum += betaBinomialPdf(i, n, alpha, beta);
        }

        return Math.min(1.0, sum);
    }

    public static double betaBinomialPvalueFast(int k, int n, double alpha, double beta) {
        return betaBinomialPvalueFast(k, n, alpha, beta, 1e-10);
    }

    public static double betaBinomialPvalueFast(int k, int n, double alpha, double beta, double pThreshold) {
        double sum = 1.0 + 0.5 * betaBinomialPdf(k, n, alpha, beta);

        for (int i = 0; i < k; i++) {
            sum -= betaBinomialPdf(i, n, alpha, beta);

            if (sum <= pThreshold) {
                return pThreshold;
            }
        }

        return sum;
    }

    private static final double SQRT2 = Math.sqrt(2.0);

    public static double normalCdf(double x, double mean, double sd) {
        final double dev = x - mean;
        if (FastMath.abs(dev) > 40 * sd) {
            return dev < 0 ? 0.0d : 1.0d;
        }
        return 0.5 * Erf.erfc(-dev / (sd * SQRT2));
    }

    public static double hypgeomCdf(int n12, int n1, int n2, int nSamples) {
        // http://journals.sagepub.com.sci-hub.cc/doi/pdf/10.2466/pms.1998.87.1.51

        int v = Math.max(0, n1 + n2 - nSamples),
                w = Math.min(n1, n2);

        assert n12 >= v && n12 <= w;

        double pPrev = 1, T = pPrev, S = pPrev;

        for (int i = v + 1; i <= w; i++) {
            pPrev *= ((double) (n1 - i + 1) * (n2 - i + 1)) / i / (nSamples - n1 - n2 + i);
            T += pPrev;

            if (i == n12) {
                S = T - 0.5 * pPrev;
            }
        }

        return S / T;
    }

    public static double hypgeomPValue(int n12, int n1, int n2, int nSamples) {
        double p = hypgeomCdf(n12, n1, n2, nSamples);

        return Math.min(p, 1 - p);
    }
}