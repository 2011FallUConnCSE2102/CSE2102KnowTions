/*
 * (C) Copyright 2005, Gregor Heinrich (gregor :: arbylon : net) (This file is
 * part of the lda-j (org.knowceans.lda.*) experimental software package.)
 */
/*
 * lda-j is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 */
/*
 * lda-j is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 */
/*
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*
 * Created on Jan 4, 2005
 */
package org.knowceans.lda;

import static java.lang.Math.exp;
import static java.lang.Math.log;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.knowceans.util.MatrixIo;

/**
 * Utility functions
 * <p>
 * lda-c reference: functions in utils.c
 * 
 * @author heinrich
 */
public class Utils {
    /*
     * given log(a) and log(b), return log(a + b)
     */
    public static double logSum(double log_a, double log_b) {
        double v;

        if (log_a < log_b) {
            v = log_b + log(1 + exp(log_a - log_b));
        } else {
            v = log_a + log(1 + exp(log_b - log_a));
        }
        return (v);
    }

    /*
     * taylor approximation of first derivative of the log gamma function
     */

    public static double digamma(double x) {
        double p;
        assert x > 0;
        x = x + 6;
        p = 1 / (x * x);
        p = (((0.004166666666667 * p - 0.003968253986254) * p + 0.008333333333333)
            * p - 0.083333333333333)
            * p;
        p = p + log(x) - 0.5 / x - 1 / (x - 1) - 1 / (x - 2) - 1 / (x - 3) - 1
            / (x - 4) - 1 / (x - 5) - 1 / (x - 6);
        return p;
    }

    static boolean GAMMAACCURATE = true;

    static final double HALFLN2PI = 9.189385332046727E-001;

    static final double HALFLN2 = 3.465735902799726e-001;

    static final double INV810 = 1.234567901234568E-003;

    public static double lgamma(double x) {
        double lnx, einvx;
        double prec;
        assert x > 0;
        lnx = log(x);
        einvx = exp(1. / x);

        if (GAMMAACCURATE) {
            prec = x * x * x;
            prec *= prec;
            prec = INV810 / prec;
            /*
             * y = x * ( log(x) - 1 + .5 * log(x * sinh(1/x) + prec) ) - .5 *
             * log(x) + .5 * log(2 * pi)
             */
            return x
                * (lnx - 1. + .5 * log(x * (einvx - 1. / einvx) / 2. + prec))
                - .5 * lnx + HALFLN2PI;
        }
        /*
         * y = x * ( 1.5 * log(x) - 1 + .5 * log(exp(1/x) - 1/exp(1/x)) - .5 *
         * log(2) ) - .5 * log(x) + .5 * log(2 * pi);
         */
        return x * (1.5 * lnx - 1. + .5 * log(einvx - 1. / einvx) - HALFLN2)
            - .5 * lnx + HALFLN2PI;
    }

    public static double fgamma(double x) {
        return Math.exp(lgamma(x));
    }

    public static long faculty(long n) {
        return (long) fgamma(n + 1);
    }

    public static long faculty(int n) {
        return (long) fgamma(n + 1);
    }

    public static double log_gamma(double x) {
        double z;
        assert x > 0;
        z = 1. / (x * x);

        x = x + 6;
        z = (((-0.000595238095238 * z + 0.000793650793651) * z - 0.002777777777778)
            * z + 0.083333333333333)
            / x;
        z = (x - 0.5) * log(x) - x + 0.918938533204673 + z - log(x - 1)
            - log(x - 2) - log(x - 3) - log(x - 4) - log(x - 5) - log(x - 6);
        return z;
    }

    static NumberFormat nf = new DecimalFormat("0.00000");

    /**
     * @param d
     * @return
     */
    public static String formatDouble(double d) {
        String x = nf.format(d);
        // String x = shadeDouble(d, 1);
        return x;

    }

    public static void main(String[] args) {
        double[][] a = { {25, 4, 1, 5}, {34, 3, 23, 55}};
        MatrixIo.saveBinaryMatrix("test.zip", a);
        double[][] b = MatrixIo.loadBinaryMatrix("test.zip");
        for (int i = 0; i < b.length; i++) {
            for (int j = 0; j < b[i].length; j++) {
                if (j > 0)
                    System.out.print(" ");
                System.out.print(b[i][j]);
            }
            System.out.println();
        }
    }
}
