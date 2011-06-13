/*
 * Created on Jun 1, 2004 To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
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

import java.util.Arrays;

import javax.media.jai.Histogram;

/**
 * Sampling methods for beta, gamma, multinomial, and Dirichlet distributions as
 * well as Dirichlet processes, using Sethurahman's stick-breaking construction
 * and Chinese restaurant process.
 * <p>
 * FIXME: markov condition in random generator, see random string?
 * 
 * @author heinrich (partly adapted from Yee Whye Teh's npbayes Matlab / C code)
 */
public class Samplers {

    private static boolean haveNextNextGaussian = false;

    private static double nextNextGaussian;

    static double drand48() {
        return Cokus.randDouble();
    }

    /**
     * uses same approach as java.util.Random()
     * 
     * @param mu
     * @param sigma
     * @return
     */
    public static double randNorm(double mu, double sigma) {
        if (haveNextNextGaussian) {
            haveNextNextGaussian = false;
            return nextNextGaussian;
        } else {
            double v1, v2, s;
            do {
                v1 = 2 * drand48() - 1; // between -1 and 1
                v2 = 2 * drand48() - 1; // between -1 and 1
                s = v1 * v1 + v2 * v2;
            } while (s >= 1 || s == 0);
            double multiplier = Math.sqrt(-2 * Math.log(s) / s);
            nextNextGaussian = v2 * multiplier;
            haveNextNextGaussian = true;
            return v1 * multiplier * sigma + mu;
        }
    }

    /**
     * GMM sampling
     * 
     * @param probs mixture responsibilities
     * @param mean mean vector
     * @param sigma stddev vector
     * @return
     */
    public static double[] randGmm(int n, double[] probs, double[] mean,
        double[] sigma) {
        double[] x = new double[n];
        int k = probs.length;
        // init random number generator

        // multinomial cdf for probs
        double[] cumprobs = new double[k];
        cumprobs[0] = probs[0];
        for (int i = 1; i < k; i++) {
            cumprobs[i] = cumprobs[i - 1] + probs[i];
        }

        for (int i = 0; i < n; i++) {

            // multinomial index sampling
            double s = drand48();
            int c = 0;
            for (c = 0; c < k; c++) {
                if (s < cumprobs[c])
                    break;
            }
            // normal component sampling
            x[i] = randNorm(mean[c], sigma[c]);
        }
        return x;
    }

    /**
     * beta as two-dimensional Dirichlet
     * 
     * @param aa
     * @param bb
     * @return
     */
    public static double randBeta(double aa, double bb) {

        double[] p = randDir(new double[] {aa, bb});
        return p[0];
    }

    /**
     * randbeta(aa, bb) Generates beta samples, one for each element in aa/bb,
     * and scale 1.
     * 
     * @param aa
     */
    public static double[] randBeta(double[] aa, double[] bb) {
        double[] beta = new double[aa.length];
        for (int i = 0; i < beta.length; i++) {
            beta[i] = randBeta(aa[i], bb[i]);
        }
        return beta;
    }

    /**
     * self-contained gamma generator (from npbayes-2.1). Multiply result with
     * scale parameter (or divide by rate parameter).
     * 
     * @param rr shape parameter
     * @return
     */
    public static double randGamma(double rr) {
        double aa, bb, cc, dd;
        double uu, vv, ww, xx, yy, zz;

        if (rr <= 0.0) {
            /* Not well defined, set to zero and skip. */
            return 0.0;
        } else if (rr == 1.0) {
            /* Exponential */
            return -Math.log(drand48());
        } else if (rr < 1.0) {
            /* Use Johnks generator */
            cc = 1.0 / rr;
            dd = 1.0 / (1.0 - rr);
            while (true) {
                xx = Math.pow(drand48(), cc);
                yy = xx + Math.pow(drand48(), dd);
                if (yy <= 1.0) {
                    assert yy != 0 && xx / yy > 0;
                    return -Math.log(drand48()) * xx / yy;
                }
            }
        } else { /* rr > 1.0 */
            /* Use bests algorithm */
            bb = rr - 1.0;
            cc = 3.0 * rr - 0.75;
            while (true) {
                uu = drand48();
                vv = drand48();
                ww = uu * (1.0 - uu);
                yy = Math.sqrt(cc / ww) * (uu - 0.5);
                xx = bb + yy;
                if (xx >= 0) {
                    zz = 64.0 * ww * ww * ww * vv * vv;
                    assert zz > 0 && bb != 0 && xx / bb > 0;
                    if ((zz <= (1.0 - 2.0 * yy * yy / xx))
                        || (Math.log(zz) <= 2.0 * (bb * Math.log(xx / bb) - yy))) {
                        return xx;
                    }
                }
            }
        }
    }

    /**
     * randgamma(aa) Generates gamma samples, one for each element in aa.
     * 
     * @param aa
     */
    public static double[] randGamma(double[] aa) {
        double[] gamma = new double[aa.length];
        for (int i = 0; i < gamma.length; i++) {
            gamma[i] = randGamma(aa[i]);
        }
        return gamma;
    }

    /**
     * Random permutation of size elements (symbols '0'.. '[size-1]'). This
     * works a bit like sampling without replacement or a factorial.
     * 
     * @param size
     * @return
     */
    public static int[] randPerm(int size) {
        int[] perm = Vectors.range(0, size - 1);
        for (int i = size - 1; i > 0; i--) {
            int k = (int) (Cokus.randDouble() * (i + 1));
            if (k != i) {
                int buf = perm[i];
                perm[i] = perm[k];
                perm[k] = buf;
            }
        }
        return perm;
    }

    /**
     * symmetric Dirichlet sample.
     * 
     * @param aa
     * @return
     */
    public static double[] randDir(double a, int dimension) {
        double[] aa = new double[dimension];
        Arrays.fill(aa, a);
        return randDir(aa);
    }

    /**
     * randdir(aa) generates one Dirichlet sample vector according to the
     * parameters alpha. ORIG: Generates Dirichlet samples, with weights given
     * in aa. The output sums to 1 along normdim, and each such sum corresponds
     * to one Dirichlet sample.
     * 
     * @param aa
     * @param normdim
     * @return
     */
    public static double[] randDir(double[] aa) {
        // function ww = randdir(aa,normdim);
        //
        // ww = randgamma(aa);
        double[] ww = randGamma(aa);
        // index.subs{normdim} = ones(1,size(aa,normdim));

        // % Dir(vecAlpha) = 1/Sum(Gamma(vecAlpha)) * Gamma(vecAlpha)
        // ww = ww./subsref(sum(ww,normdim),index);
        double sum = 0;
        for (int i = 0; i < ww.length; i++) {
            sum += ww[i];
        }
        for (int i = 0; i < ww.length; i++) {
            ww[i] /= sum;
        }
        return ww;
    }

    /**
     * Generate as many Dirichlet column samples as there are columns (direction =
     * 1; randdir(A, 1)) or row samples as there are rows (direction = 2,
     * randdir(A, 2)) in aa (aa[][]), taking the respective parameters.
     * 
     * @param aa
     * @param direction -- 2 is more efficient (row-major Java matrix structure)
     * @return
     */
    public static double[][] randDir(double[][] aa, int direction) {
        double[][] ww = null;
        if (direction == 1) {
            ww = new double[aa.length][aa[0].length];
            double[] dirsmp = new double[aa.length];
            for (int i = 0; i < ww.length; i++) {
                for (int j = 0; j < dirsmp.length; j++) {
                    dirsmp[i] = aa[j][i];
                }
                dirsmp = randDir(aa[i]);
                for (int j = 0; j < dirsmp.length; j++) {
                    ww[j][i] = dirsmp[j];
                }
            }
        } else {
            ww = new double[aa.length][aa[0].length];
            for (int i = 0; i < ww.length; i++) {
                ww[i] = randDir(aa[i]);
            }
        }
        return ww;
    }

    /**
     * Generate n Dirichlet samples taking parameters aa.
     * 
     * @param aa
     * @return
     */
    public static double[][] randDir(double[] aa, int repetitions) {
        double[][] ww = new double[repetitions][aa.length];
        for (int i = 0; i < repetitions; i++) {
            ww[i] = randDir(aa);
        }
        return ww;
    }

    /**
     * Multiply sample a multinomial distribution and return a vector with
     * category frequencies.
     * 
     * @param pp
     * @param repetitions
     * @return vector of frequencies of the categories
     */
    public static int[] randMultFreqs(double[] pp, int repetitions) {
        int[] freqs = new int[pp.length];
        for (int i = 0; i < freqs.length; i++) {
            freqs[i] = 0;
        }

        for (int i = 0; i < repetitions; i++) {
            freqs[randMult(pp)]++;
        }

        return freqs;
    }

    /**
     * Multiply sample a multinomial distribution and return a vector with all
     * samples.
     * 
     * @param pp
     * @param repetitions
     * @return vector of all samples.
     */
    public static int[] randMult(double[] pp, int repetitions) {
        int[] samples = new int[repetitions];

        for (int i = 0; i < repetitions; i++) {
            samples[i] = randMult(pp);
        }
        return samples;

    }

    public static int randMultSimple(final double[] pp) {

        // pp = cumsum(pp);
        // ii = 1+sum(rand*pp(end)>pp);
        int i;
        double[] cumPp = new double[pp.length];

        System.arraycopy(pp, 0, cumPp, 0, pp.length);

        for (i = 1; i < pp.length; i++) {
            cumPp[i] += cumPp[i - 1];

        }
        // this automatically normalises.
        double randNum = Cokus.randDouble() * cumPp[i - 1];

        // TODO: use binarySearch().
        for (i = 0; i < cumPp.length; i++) {
            if (cumPp[i] > randNum) {
                break;
            }
        }

        return i;
    }

    public static void testMult() {
        int N = 100000;

        int[][] a = new int[2][N];

        // double[] p = new double[] {.3, .4, .1, .2};
        double[] p = randDir(.1, 10000);
        System.out.println(Vectors.print(p));

        Cokus.seed(4357);
        long t0 = System.currentTimeMillis();
        for (int i = 0; i < N; i++) {
            a[0][i] = randMult(p);
        }
        long t1 = System.currentTimeMillis();
        System.out.println();
        Cokus.seed(4357);

        long t2 = System.currentTimeMillis();
        for (int i = 0; i < N; i++) {
            a[1][i] = randMult(p);
        }
        long t3 = System.currentTimeMillis();

        // for (int i = 0; i < a[0].length; i++) {
        // System.out.println(a[0][i] + " " + a[1][i]);
        // if (a[0][i] != a[1][i]) {
        // System.err.println("a0 != a1");
        // }
        // }
        System.out.println("linear search: " + (t1 - t0) + " binary search: "
            + (t3 - t2));
    }

    /**
     * Creates one multinomial sample given the parameter vector pp. Each
     * category is named after the index (0-based!) of the respective element of
     * pp; Sometimes called categorical distribution (e.g., in BUGS). This
     * version uses a binary search algorithm and does not require
     * normalisation.
     */
    public static int randMult(final double[] pp) {

        // pp = cumsum(pp);
        // ii = 1+sum(rand*pp(end)>pp);
        int i;
        double[] cumPp = new double[pp.length];

        System.arraycopy(pp, 0, cumPp, 0, pp.length);

        for (i = 1; i < pp.length; i++) {
            cumPp[i] += cumPp[i - 1];

        }
        // this automatically normalises.
        double randNum = Cokus.randDouble() * cumPp[i - 1];

        // TODO: use insertion point formula in Array.binarySearch()
        i = binarySearch(cumPp, randNum);

        return i;
    }

    /**
     * perform a binary search and return the first index i at which a[i] >= p.
     * Adapted from java.util.Arrays.binarySearch.
     * 
     * @param a
     * @param p
     * @return
     */
    public static int binarySearch(double[] a, double p) {
        if (p < a[0]) {
            return 0;
        }
        int low = 0;
        int high = a.length - 1;
        while (low <= high) {
            int mid = (low + high) >> 1;
            double midVal = a[mid];

            int cmp;
            if (midVal < p) {
                low = mid + 1;
            } else if (midVal > p) {
                if (a[mid - 1] < p)
                    return mid;
                high = mid - 1;
            } else {
                return mid;
            }
        }
        return a.length; // key not found.
    }

    /**
     * randconparam(alpha,numdata,numclass,aa,bb,numiter) Generates a sample
     * from a concentration parameter of a HDP with gamma(aa,bb) prior, and
     * number of classes and data items given in numdata, numclass (has to be
     * row vectors). Uses auxiliary variable method, for numiter iterations.
     * <p>
     * Modification of Escobar and West. Works for multiple groups of data.
     * numdata, numclass are row vectors, one element per group. from
     * npbayes-2.1
     * 
     * @param alpha alpha
     * @param numgroup number of components ??
     * @param numdata number of data items per class
     * @param numtable number of per DP
     * @param alphaa hyperparameter (gamma shape)
     * @param alphab hyperparameter (gamma scale)
     * @param numiter number of iterations
     * @return
     */
    public static double randConParam(double alpha, int numgroup,
        int[] numdata, int[] numtable, double alphaa, double alphab, int numiter) {
        int iter, jj, nd, zz;
        double aa, bb, eta;

        // Teh: Escobar and West's method for single Gamma.

        for (iter = 0; iter < numiter; iter++) {
            aa = alphaa;
            bb = alphab;
            for (jj = 0; jj < numgroup; jj++) {
                nd = numdata[jj];
                eta = randBeta(alpha + 1.0, nd);
                // zz = (drand48() * (alpha + nd) < nd);
                zz = (drand48() * (alpha + nd) < nd) ? 1 : 0;

                aa += numtable[jj] - (zz);
                bb -= Math.log(eta);
            }
            alpha = randGamma(aa) / bb;
        }
        return alpha;
    }

    /**
     * Sample the Dirichlet process concetration parameter given the topic and
     * data counts and gamma hyperparameters alphaa and alphab.
     * 
     * @param alpha
     * @param numdata
     * @param numtopic
     * @param alphaa
     * @param alphab
     * @param numiter
     * @return
     */
    public static double randConParam(double alpha, int numdata, int numtopic,
        double alphaa, double alphab, int numiter) {
        int iter, jj, nd, zz;
        double aa, bb, eta;

        // Escobar and West's method
        double alphaAvg = 0.;

        for (iter = 0; iter < numiter; iter++) {
            aa = alphaa;
            bb = alphab;

            // e+w (14)
            eta = randBeta(alpha + 1.0, numdata);

            // e+w (13)
            double pi = 1 / (numdata * (alphab - Math.log(eta))
                / (alphaa + numtopic - 1) + 1);
            // choose between the two gamma components
            zz = (drand48() > pi) ? 1 : 0;

            // sample gamma
            aa += numtopic - zz;
            bb -= Math.log(eta);
            // * or / ? (scale or rate?)
            alpha = randGamma(aa) / bb;
            alphaAvg += alpha;
        }
        return alphaAvg / numiter;
    }

    // /**
    // * randconparam(alpha,numdata,numclass,aa,bb,numiter) Generates a sample
    // * from a concentration parameter with gamma(aa,bb) prior, and number of
    // * classes and data items given in numdata, numclass (has to be row
    // * vectors). Uses auxiliary variable method, for numiter iterations.
    // * <p>
    // * Escobar and West's method for single Gamma.
    // */
    // public static double randConParam1(double alpha, int numdata, int
    // numclass,
    // double aa, double bb, int numiter) {
    // // function alpha = randconparam(alpha,numdata,numclass,aa,bb,numiter);
    // //
    // // % Escobar and West's method for single Gamma.
    // //
    // // if nargin == 5
    // // numiter = 1;
    // // end
    // //
    // // for ii = 1:numiter
    // // xx = randbeta(alpha+1,numdata);
    // //
    // // weights = [(aa+numclass-1)/(bb-log(xx)) numdata];
    // // zz = randmult(weights/sum(weights));
    // //
    // // if zz==1
    // // alpha = randgamma(aa+numclass,1./(bb-log(xx)));
    // // else
    // // alpha = randgamma(aa+numclass-1,1./(bb-log(xx)));
    // // end
    // // end
    // return 0;
    // }
    //
    // /**
    // * randconparam(alpha,numdata,numclass,aa,bb,numiter) Generates a sample
    // * from a concentration parameter with gamma(aa,bb) prior, and number of
    // * classes and data items given in numdata, numclass (has to be row
    // * vectors). Uses auxiliary variable method, for numiter iterations.
    // * <p>
    // * less stable for small alpha's as xx gets very close to zero.
    // */
    // public static double randConParam2(double alpha, int numdata, int
    // numclass,
    // double aa, double bb, int numiter) {
    // // function alpha = randconparam2(alpha,numdata,numclass,aa,bb);
    // //
    // // % less stable for small alpha's as xx gets very close to zero.
    // //
    // // xx = randbeta(alpha,numdata);
    // //
    // // alpha = randgamma(aa+numclass,1./(bb-log(xx)));
    // //
    // return 0;
    // }
    //
    // /**
    // * randconparam(alpha,numdata,numclass,aa,bb,numiter) Generates a sample
    // * from a concentration parameter with gamma(aa,bb) prior, and number of
    // * classes and data items given in numdata, numclass (has to be row
    // * vectors). Uses auxiliary variable method, for numiter iterations.
    // * <p>
    // * Modification of Escobar and West. Works for multiple groups of data.
    // * numdata, numclass are row vectors, one element per group.
    // */
    // public static double randConParam3(double alpha, int numdata, int
    // numclass,
    // double aa, double bb, int numiter) {
    // // function alpha = randconparam3(alpha,numdata,numclass,aa,bb,numiter);
    // //
    // // % Modification of Escobar and West. Works for multiple groups of
    // // data.
    // // % numdata, numclass are row vectors, one element per group.
    // //
    // // if nargin == 5
    // // numiter = 1;
    // // end
    // //
    // // totalclass = sum(numclass);
    // // num = length(numdata);
    // //
    // // for ii = 1:numiter
    // // xx = randbeta((alpha+1)*ones(1,num),numdata);
    // //
    // // zz = rand(1,num).*(alpha+numdata)<numdata;
    // //
    // // gammaa = aa + totalclass - sum(zz);
    // // gammab = bb - sum(log(xx));
    // // alpha = randgamma(gammaa)./gammab;
    // //
    // // end
    // return 0;
    // }

    public static CrpData randCrp(double alpha, int numdata) {
        return randCrp(new double[] {alpha}, numdata);
    }

    /**
     * [cc numclass] = randcrp(alpha,numdata) Generates a partition of numdata
     * items with concentration parameter alpha, which can be an array, in which
     * case the Chinese restaurant process has "two new tables to chose for each
     * new customer". cc is sequence of indicator variables denoting which class
     * each data item is in ("on which table each customer sits"), and numclass
     * is the generated number of classes.
     * 
     * @param alpha
     * @param numdata
     * @return
     */
    public static CrpData randCrp(double[] alpha, int numdata) {
        // function [cc, numclass] = randcrp(alpha, numdata);
        // % generates a CRP partition of numdata items, with concentration
        // parameter
        // % alpha.
        //
        // cc = zeros(1,numdata);
        // weights = alpha;
        // ? vector

        double[] weights = Vectors.copy(alpha);
        // numclass = 0;
        CrpData crp = new CrpData(numdata);
        //
        // for ii = 1:numdata
        // cc(ii) = randmult(weights);
        // if cc(ii) > numclass
        // weights = [weights(1:numclass) 1 alpha];
        // numclass = cc(ii);
        // else
        // weights(cc(ii)) = weights(cc(ii)) + 1;
        // end
        // end

        for (int ii = 0; ii < numdata; ii++) {
            crp.cc[ii] = randMult(weights);
            if (crp.cc[ii] > crp.numclass) {
                // has the multinomial weights (Vector better)
                // add one component weight (1)
                weights = Vectors.concat(Vectors.subVector(weights, 0,
                    crp.numclass - 1), new double[] {1}, alpha);
                crp.numclass = crp.cc[ii];
            } else {
                weights[crp.cc[ii]]++;
            }
        }
        return crp;
    }

    /**
     * data structure for a Chinese restaurant process CrpData
     * 
     * @author heinrich
     */
    public static class CrpData {
        public int[] cc;

        public int numclass;

        public CrpData(int numdata) {
            cc = Vectors.ones(numdata, 0);
            numclass = 0;
        }
    }

    /**
     * randnumtable(weights,maxtable) For each entry in weights and maxtables,
     * generates the number of tables given concentration parameter (weights)
     * and number of data items (maxtable). From npbayes-2.1. enumClass seems to
     * be the expected value of randNumTable
     * 
     * @param weights
     * @param maxtable
     * @return
     */
    public static int randNumTable(double alpha, int numdata) {
        int ii, numtable;

        if (numdata == 0) {
            return 0;
        } else {
            numtable = 1;
            for (ii = 1; ii < numdata; ii++) {
                if (drand48() < alpha / (ii + alpha))
                    numtable++;
            }
            return numtable;
        }
    }

    /**
     * randstick(alpha,numclass) Generates stick-breaking weights with
     * concentration parameter for numclass "sticks". XXX: untested
     * 
     * @param alpha
     * @param numclass
     * @return
     */
    public static double[] randStick(double[] alpha, int numclass) {
        // function beta = randstick(alpha,numclass);
        //
        // one = ones(1,numclass);
        // zz = randbeta(one, alpha*one);
        // beta = zz .* cumprod([1 1-zz(1:numclass-1)]);

        double[] beta = new double[numclass];
        double[] one = Vectors.ones(numclass, 1.0);
        double[] zz = new double[numclass];
        for (int i = 0; i < zz.length; i++) {
            zz[i] = randBeta(1, alpha[i]);
        }
        beta[0] = 1;
        for (int i = 1; i < numclass; i++) {
            beta[i] = (1 - zz[i - 1]) * beta[i - 1];
        }
        for (int i = 0; i < numclass; i++) {
            beta[i] *= zz[i];
        }
        return beta;
    }

    /**
     * enumclass(alpha,numdata) The expected number of tables in a CRP with
     * concentration parameter alpha and numdata items.
     * 
     * @param alpha
     * @param numdata
     * @return
     */
    public static double enumClass(double alpha, int numdata) {
        // function numclass = enumclass(alpha,numdata);
        //
        // numclass = alpha*sum(1./(alpha-1+(1:numdata)));
        double numclass = 0;
        for (int ii = 1; ii <= numdata; ii++) {
            numclass += 1 / ((alpha - 1) + ii);
        }
        return numclass * alpha;
    }

    /**
     * create a random string of length alphanumeric characters.
     * 
     * @param length of output
     * @param alphabet alphabet to be used or null
     * @return
     */
    public static String randString(int length, byte[] alphabet) {
        if (alphabet == null)
            alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"
                .getBytes();
        byte[] pass = new byte[length];

        for (int k = 0; k < length; k++) {
            int i = (int) Math.floor(Cokus.randDouble() * alphabet.length);
            pass[k] = alphabet[i];
        }
        return new String(pass);
    }

    /**
     * meanlik(lik) Computes estimated likelihood from individual samples.
     * Basically does a harmonic mean of lik in 3rd dimension, followed by
     * normal mean in 2nd.
     * 
     * @param lik
     * @return
     */
    public static double meanLik(double lik) {
        // function lik = meanlik(lik);
        //
        // lik = logmeanexp(-logmeanexp(-lik,3),2);

        return 0;
    }

    /**
     * TODO: MAXSTIRLING should be made variable
     */
    private static int MAXSTIRLING = 40;

    /**
     * maximum stirling number in allss
     */
    static int maxnn = 1;

    /**
     * contains all stirling number iteratively calculated so far
     */
    static double[][] allss = new double[MAXSTIRLING][];

    /**
     * 
     */
    static double[] logmaxss = new double[MAXSTIRLING];

    static double lmss = 0;

    /**
     * [ss lmss] = stirling(nn) Gives unsigned Stirling numbers of the first
     * kind s(nn,*) in ss. ss(i) = s(nn,i-1). ss is normalized so that maximum
     * value is 1, and the log of normalization is given in lmss (static
     * variable).
     * 
     * @param nn
     * @return
     */
    public static double[] stirling(int nn) {
        // function [ss, lmss] = stirling(nn);
        //
        // % returns the unsigned stirling numbers of the first kind.
        // % ss(nn,tt) for tt=1:nn
        //
        // persistent maxnn allss logmaxss
        //
        // if isempty(maxnn)
        // maxnn = 1;
        // allss = {1};
        // logmaxss = 0;
        // end
        //
        // if nn > maxnn
        // allss{nn} = [];
        // logmaxss(nn) = 0;
        // for mm=maxnn+1:nn
        // allss{mm} = [allss{mm-1}*(mm-1) 0] + [0 allss{mm-1}];
        // mss = max(allss{mm});
        // allss{mm} = allss{mm}/mss;
        // logmaxss(mm) = logmaxss(mm-1) + log(mss);
        // end
        // maxnn = nn;
        // end
        //
        // ss = allss{nn};
        // lmss = logmaxss(nn);

        // nn is a mathematical argument; index into allss is nn - 1
        if (allss[0] == null) {
            allss[0] = new double[1];
            allss[0][0] = 1;
            logmaxss[0] = 0;
        }

        if (nn > maxnn) {
            for (int mm = maxnn; mm < nn; mm++) {
                int len = allss[mm - 1].length + 1;
                allss[mm] = new double[len];
                for (int xx = 0; xx < len; xx++) {
                    // allss{mm} = [allss{mm-1}*(mm-1) 0] + ...
                    allss[mm][xx] += (xx < len - 1) ? allss[mm - 1][xx] * mm
                        : 0;
                    // [0 allss{mm-1}];
                    allss[mm][xx] += (xx == 0) ? 0 : allss[mm - 1][xx - 1];
                }
                double mss = Vectors.max(allss[mm]);
                Vectors.mult(allss[mm], 1 / mss);
                logmaxss[mm] = logmaxss[mm - 1] + Math.log(mss);
            }
            maxnn = nn;
        }
        lmss = logmaxss[nn - 1];
        //
        return allss[nn - 1];

    }

    /**
     * @param numclass
     * @return
     */
    public static int randUniform(int numvalue) {
        return (int) Math.floor(drand48() * (double) numvalue);
    }

    /* ********** colt-based stuff for reference *************** */

    // public static RandomEngine rand = new MersenneTwister(new Date());
    //
    // /**
    // * randbeta(aa,bb) Generates one beta sample, with weights given by aa and
    // * bb. A beta sample is the same as a 2-dimensional Dirichlet sample.
    // */
    // public static double randbeta(double aa, double bb) {
    //
    // Beta g = new Beta(aa, bb, rand);
    // return g.nextDouble();
    // }
    //
    // /**
    // * randgamma(aa) Generates one gamma sample with scale 1.
    // *
    // * @param aa
    // */
    // public static double randgamma(double aa) {
    //
    // Gamma g = new Gamma(aa, 1.0, rand);
    // return g.nextDouble();
    // }
    //
}
