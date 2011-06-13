/*
 * Created on Jul 18, 2005
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

/**
 * Densities calculates for different density functions the likelihood of a data
 * item given the parameters
 * 
 * @author heinrich
 */
public class Densities {

    /**
     * Normal likelihood
     * 
     * @param x
     * @param mu
     * @param sigma
     * @return
     */
    public static double pdfNorm(double x, double mu, double sigma) {
        double p = 1 / (Math.sqrt(2 * Math.PI) * sigma)
            / Math.exp((x - mu) * (x - mu) / (2 * sigma * sigma));
        return p;
    }

    /**
     * GMM likelihood
     * 
     * @param x
     * @param k
     * @param probs
     * @param mean
     * @param sigma
     * @return
     */
    public static double pdfGmm(double x, int k, double[] probs, double[] mean,
        double[] sigma) {
        double p = 0;
        for (int i = 0; i < k; i++) {
            p += probs[i] * pdfNorm(x, mean[i], sigma[i]);
        }
        return p;
    }

    /**
     * gamma likelihood p(x|a,b) = x^(a-1) * e^(-x/b) / (b^a * Gamma(a))
     * 
     * @param x value
     * @param a (shape?)
     * @param b (scale?)
     * @return
     */
    public static double pdfGamma(double x, double a, double b) {
        return Math.pow(x, a - 1) * Math.exp(-x / b)
            / (Math.pow(b, a) * Math.exp(Gamma.lgamma(a)));
    }

    /**
     * beta likelihood
     * 
     * @param x data item
     * @param a pseudo counts for success
     * @param b pseudo counts for failure
     * @return
     */
    public static double pdfBeta(double x, double a, double b) {
        return pdfDirichlet(new double[] {x, 1 - x}, new double[] {a, b});
    }

    /**
     * Dirichlet likelihood using logarithmic calculation
     * <p>
     * Dir(xx|alpha) =
     * <p>
     * Gamma(sum_i alpha[i])/(prod_i Gamma(alpha[i])) prod_i xx[i]^(alpha[i]-1)
     * 
     * @param xx multivariate convex data item (sum=1)
     * @param alpha Dirichlet parameter vector
     * @return
     */
    public static double pdfDirichlet(double[] xx, double[] alpha) {
        double sumAlpha = 0.;
        double sumLgammaAlpha = 0.;
        double logLik = 0.;
        for (int i = 0; i < xx.length; i++) {
            sumAlpha += alpha[i];
            sumLgammaAlpha += Gamma.lgamma(alpha[i]);
            logLik += Math.log(xx[i]) * (alpha[i] - 1);
        }
        return Math.exp(Gamma.lgamma(sumAlpha) - sumLgammaAlpha + logLik);
    }

    /**
     * Symmetric Dirichlet likelihood:
     * <p>
     * Dir(xx|alpha) = Gamma(k * alpha)/Gamma(alpha)^k prod_i xx[i]^(alpha - 1)
     * 
     * @param xx multivariate convex data item (sum=1)
     * @param alpha symmetric parameter
     * @return
     */
    public static double pdfDirichlet(double[] xx, double alpha) {
        double logCoeff = Gamma.lgamma(alpha * xx.length) - Gamma.lgamma(alpha)
            * xx.length;
        double logLik = 0.;
        for (int i = 0; i < xx.length; i++) {
            logLik += Math.log(xx[i]);
        }
        logLik *= alpha - 1;
        return Math.exp(logCoeff + logLik);
    }

    /**
     * Mult(nn|pp) using logarithmic multinomial coefficient
     * 
     * @param nn counts for each category
     * @param pp convex probability vector for categories
     * @return
     */
    public static double pdfMultinomial(int[] nn, double[] pp) {
        int N = 0;
        double logCoeff = 0.;
        double logLik = 0.;
        for (int i = 0; i < nn.length; i++) {
            N += nn[i];
            logCoeff -= Gamma.lgamma(nn[i] + 1);
            logLik += Math.log(pp[i]) * nn[i];
        }
        logCoeff += Gamma.lgamma(N + 1);
        return Math.exp(logCoeff + logLik);
    }

    /**
     * Binom(n | N, p) using linear binomial coefficient
     * 
     * @param n
     * @param N
     * @param p
     * @return
     */
    public static double pdfBinomial(int n, int N, double p) {
        long binom = Gamma.faculty(N)
            / (Gamma.faculty(N - n) * Gamma.faculty(n));
        double lik = binom * Math.pow(p, N) * Math.pow(1 - p, N - n);
        return lik;
    }

    public static void main(String[] args) {
        double[] alpha = new double[] {2, 1};
        for (float i = 0; i < 1; i += .01) {
            System.out.println(i + "   "
                + pdfDirichlet(new double[] {i, 1 - i}, alpha));
        }

    }
}
