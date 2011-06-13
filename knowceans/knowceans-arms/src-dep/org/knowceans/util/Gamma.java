/*
 * Created on Aug 1, 2005
 */
/*
 * (C) Copyright 2005, Gregor Heinrich (gregor :: arbylon : net) (This file is
 * part of the knowceans-tools (org.knowceans.util.*) experimental software package.)
 */
/*
 * knowceans-tools is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 */
/*
 * knowceans-tools is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 */
/*
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.knowceans.util;

import static java.lang.Math.log;

/**
 * Gamma represents the Gamma function and its derivatives
 * 
 * @author heinrich
 */
public class Gamma {

    /**
     * truncated Taylor series of log Gamma(x). From lda-c
     * 
     * @param x
     * @return
     */
    public static double lgamma(double x) {
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

    /**
     * gamma function
     * 
     * @param x
     * @return
     */
    public static double fgamma(double x) {
        return Math.exp(lgamma(x));
    }

    /**
     * faculty of an integer.
     * 
     * @param n
     * @return
     */
    public static int faculty(int n) {
        return (int) Math.exp(lgamma(n + 1));
    }

    /**
     * "Dirichlet delta function" is the partition function of the Dirichlet
     * distribution and the k-dimensional generalisation of the beta function.
     * fdelta(a) = prod_k fgamma(a_k) / fgamma( sum_k a_k ) = int_(sum x = 1)
     * prod_k x_k^(a_k-1) dx. See G. Heinrich: Parameter estimation for text
     * analysis (http://www.arbylon.net/publications/text-est_iv.pdf)
     * 
     * @param x
     * @return
     */
    public static double fdelta(double[] x) {
        double lognum = 1;
        double den = 0;
        for (int i = 0; i < x.length; i++) {
            lognum += lgamma(x[i]);
            den += x[i];
        }
        return Math.exp(lognum - lgamma(den));
    }

    /**
     * truncated Taylor series of Psi(x) = d/dx Gamma(x). From lda-c
     * 
     * @param x
     * @return
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

    /**
     * truncated Taylor series of d/dx Psi(x) = d^2/dx^2 Gamma(x). From lda-c
     * 
     * @param x
     * @return
     */
    public static double trigamma(double x) {
        double p;
        int i;

        x = x + 6;
        p = 1 / (x * x);
        p = (((((0.075757575757576 * p - 0.033333333333333) * p + 0.0238095238095238)
            * p - 0.033333333333333)
            * p + 0.166666666666667)
            * p + 1)
            / x + 0.5 * p;
        for (i = 0; i < 6; i++) {
            x = x - 1;
            p = 1 / (x * x) + p;
        }
        return (p);
    }
}
