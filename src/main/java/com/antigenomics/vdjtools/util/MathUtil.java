/*
 * Copyright 2013-2014 Mikhail Shugay (mikhail.shugay@gmail.com)
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
 * Last modified on 10.11.2014 by mikesh
 */

package com.antigenomics.vdjtools.util;

import java.util.Random;

public class MathUtil {
    public double distanceCorrelation(double[][] x, double[][] y) {
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

    private static double euclDist(double[] xx, double[] yy) {
        int n = xx.length;
        if (n != yy.length)
            throw new IndexOutOfBoundsException("x and y should be of same length");

        double dist = 0;

        for (int i = 0; i < n; i++) {
            dist += xx[i] * yy[i];
        }

        return Math.sqrt(dist);
    }

    public static <T> void shuffle(T[] arr) {
        Random rnd = new Random();
        for (int i = arr.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            T tmp = arr[index];
            arr[index] = arr[i];
            arr[i] = tmp;
        }
    }

    public static double JSD(double[] pArr, double[] qArr) throws Exception {
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
