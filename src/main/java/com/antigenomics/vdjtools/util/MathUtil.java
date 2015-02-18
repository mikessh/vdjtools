/*
 * Copyright 2013-2015 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Last modified on 30.1.2015 by mikesh
 */

package com.antigenomics.vdjtools.util;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

import java.util.Random;

/**
 * Fast Java implementation for mathematical routines used in VDJtools
 */
public class MathUtil {
    private static final double HALF_LOG_2_PI = 0.5 * FastMath.log(2.0 * FastMath.PI);
    private static final long[] FACTORIALS = new long[]{
            1l, 1l, 2l,
            6l, 24l, 120l,
            720l, 5040l, 40320l,
            362880l, 3628800l, 39916800l,
            479001600l, 6227020800l, 87178291200l,
            1307674368000l, 20922789888000l, 355687428096000l,
            6402373705728000l, 121645100408832000l, 2432902008176640000l};

    /**
     * Computes the logarithm of n!/(n-k)!. Implemented based on apache math 3
     */
    public static double logFactorialRatio(final int n, final int k) {
        if (k > n)
            throw new IllegalArgumentException("k should be less or equal to n");

        return fastLogFactorial(n) - fastLogFactorial(n - k);
    }

    private static double fastLogFactorial(final long x) {
        return x < FACTORIALS.length ? FastMath.log(FACTORIALS[(int) x]) : fastLogGamma(x);
    }

    /**
     * Internal. Implemented based on apache math 3 
     * @param x
     * @return
     */
    private static double fastLogGamma(final double x) {
        double sum = Gamma.lanczos(x);
        double tmp = x + Gamma.LANCZOS_G + .5;
        return ((x + .5) * FastMath.log(tmp)) - tmp +
                HALF_LOG_2_PI + FastMath.log(sum / x);
    }

    public static double distanceCorrelation(final double[][] x, final double[][] y) {
        int n = x.length;
        if (n != y.length)
            throw new IndexOutOfBoundsException("x and y should be of same length");

        double[][] a = new double[n][n], b = new double[n][n];
        double[] a1m = new double[n], a2m = new double[n], b1m = new double[n], b2m = new double[n];
        double am = 0, bm = 0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    double aa = euclDist(x[i], x[j]), bb = euclDist(y[i], y[j]);

                    a[i][j] = aa;
                    b[i][j] = bb;

                    a1m[i] += aa;
                    a2m[j] += aa;

                    b1m[i] += bb;
                    b2m[j] += bb;

                    am += aa;
                    bm += bb;
                }
            }
        }

        double[][] A = new double[n][n], B = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = a[i][j] - (a1m[i] - a2m[j] - am / n) / n;
                B[i][j] = b[i][j] - (b1m[i] - b2m[j] - bm / n) / n;
            }
        }

        double dCovXY = 0, dVarX = 0, dVarY = 0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double AA = A[i][j], BB = B[i][j];
                dCovXY += AA * BB;
                dVarX += AA * AA;
                dVarY += BB * BB;
            }
        }

        return Math.sqrt(dCovXY) / Math.sqrt(Math.sqrt(dVarX * dVarY));
    }

    private static double euclDist(final double[] xx, final double[] yy) {
        int n = xx.length;
        if (n != yy.length)
            throw new IndexOutOfBoundsException("x and y should be of same length");

        double dist = 0;

        for (int i = 0; i < n; i++) {
            dist += xx[i] * yy[i];
        }

        return Math.sqrt(dist);
    }

    public static <T> void shuffle(final T[] arr) {
        Random rnd = new Random();
        for (int i = arr.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            T tmp = arr[index];
            arr[index] = arr[i];
            arr[i] = tmp;
        }
    }

    public static double JSD(final double[] pArr, final double[] qArr) throws Exception {
        int n = pArr.length;

        if (n != qArr.length)
            throw new Exception("Input histograms must be of same length");

        double pSum = 0, qSum = 0;

        for (int i = 0; i < n; i++) {
            pSum += pArr[i];
            qSum += qArr[i];
        }

        double jsd = 0;
        for (int i = 0; i < n; i++) {
            double p = pArr[i] / pSum, q = qArr[i] / qSum,
                    m = (p + q) / 2.0;
            jsd += (p > 0 ? (Math.log(p / m) * p) : 0d) + (q > 0 ? (Math.log(q / m) * q) : 0d);
        }

        return jsd / 2.0 / Math.log(2.0);
    }
}
