package org.knowceans.sandbox;

import org.knowceans.util.Vectors;

/**
 * FirFilter implements a simple FIR filter, exploiting the symmetries of the
 * complex spectrum of real signals.
 * <p>
 * Conventions: suffix z specifies a complex quantity or its size, they are
 * stored with interleaved real and imaginary parts.
 * 
 * @author gregor heinrich
 * @date dec. 2007, updated jan. 2010
 */
public class FirFilter {

    final static double TWOPI = 6.28318530717959;

    /**
     * create FFT coefficients for size N
     * 
     * @param Wz N/2 complex values
     * @param Nz FFT size
     */
    public static void ffw(double[] Wz, int Nz) {
        int i;
        double arg;
        for (i = 0; i < Nz / 2; i++) {
            arg = i * TWOPI / Nz;
            Wz[2 * i] = Math.cos(arg);
            Wz[2 * i + 1] = -Math.sin(arg);
        }
    }

    /**
     * forward FFT
     * 
     * @param Xz data (complex interleaved)
     * @param Wz coefficients (complex interleaved)
     * @param Nz FFT size (complex)
     */
    public static void fft(double[] Xz, double[] Wz, int Nz) {

        // permute input
        perm(Xz, Nz);

        // perform inverse transform
        fft1d(Xz, Wz, 1, Nz);
    }

    /**
     * inverse FFT
     * 
     * @param Xz data (complex interleaved)
     * @param Wz coefficients (complex interleaved)
     * @param Nz FFT size (complex)
     */
    public static void ifft(double[] Xz, double[] Wz, int Nz) {

        // permute input
        perm(Xz, Nz);

        // perform inverse transform
        fft1d(Xz, Wz, -1, Nz);
    }

    /**
     * forward real-to-complex FFT with half length
     * 
     * @param X data [N] (real)
     * @param Yz output [N/2 + 1] (complex)
     * @param Wz [N/2] coefficients (complex interleaved, of full-size FFT)
     * @param N FFT size (real)
     */
    public static void fftr(double[] X, double[] Wz, int N) {
        //fftr1d(X, Wz, 1, N);
        //fftr(X, X, Wz, N);
    }

    /**
     * inverse complex-to-real FFT with half length
     * 
     * @param Xz data (complex interleaved)
     * @param Wz coefficients (complex interleaved)
     * @param Nz FFT size (complex)
     */
    public static void ifftr(double[] Xz, double[] Wz, int Nz) {
        // perform inverse transform
        //fftr1d(Xz, Wz, -1, Nz);
    }

    /**
     * Frequency domain filtering of two real signals with two different real
     * FIR filters. Filtering of a complex signal is possible by using H1z ==
     * H2z.
     * 
     * @param Wz [N/2] twiddling factors (complex interleaved)
     * @param Xz [N] complex input vector (pair = real, imag) (complex
     *        interleaved)
     * @param tempz [N] temporary buffer
     * @param Yz [N] complex output vector (complex interleaved)
     * @param H1z [N/2 + 1] DFT of IR for filter 1 (complex interleaved)
     * @param H2z --"-- for filter 2 (complex interleaved)
     * @param Nz length of input and output vectors (complex)
     * @return
     */
    public static long dualFir(double[] Wz, double[] Xz, double[] tempz,
        double[] Yz, double[] H1z, double[] H2z, int Nz) {

        int k, l;
        int Ny = Nz / 2;
        int Nyr = 2 * Ny;
        double x;

        // use standard fft function
        System.arraycopy(Xz, 0, tempz, 0, 2 * Nz);
        fft(tempz, Wz, Nz);

        // need to deinterleave
        if (H1z != H2z) {

            // DC and Nyquist frequency
            Yz[1] = tempz[1] * H2z[0];
            Yz[0] = tempz[0] * H1z[0];
            Yz[Nyr + 1] = tempz[Nyr + 1] * H2z[Nyr];
            Yz[Nyr] = tempz[Nyr] * H1z[Nyr];

            // symmetric FFT approach
            for (k = 1; k < Ny; k++) {
                int kr = 2 * k;
                int ki = kr + 1;
                l = Nz - k;
                int lr = 2 * l;
                int li = lr + 1;
                // Y[k] = S2[k] = -j*(Z[k] - Z'[l]) / 2
                //      = (imZ[k] + imZ[l] + j*(reZ[l] - reZ[k])) / 2;
                Yz[kr] = (tempz[ki] + tempz[li]) * 0.5;
                Yz[ki] = (tempz[lr] - tempz[kr]) * 0.5;
                // T[k] = S1[k] = (Z[k] + Z'[l]) / 2
                //      = (reZ[k] + reZ[l] + j*(imZ[k] - imZ[l])) / 2
                tempz[kr] = (tempz[kr] + tempz[lr]) * 0.5;
                tempz[ki] = (tempz[ki] - tempz[li]) * 0.5;

                // Y[l] = reY[k] - j*imY[k];
                Yz[lr] = Yz[kr];
                Yz[li] = -Yz[ki];
                // T[l] = reT[k] - j*imT[k];
                tempz[lr] = tempz[kr];
                tempz[li] = -tempz[ki];

                // freq domain filter, sequence reconstruction
                double[] c = cmult(tempz[kr], tempz[ki], H1z[kr], H1z[ki]);
                tempz[kr] = c[0];
                tempz[ki] = c[1];
                c = cmult(Yz[kr], Yz[ki], H2z[kr], H2z[ki]);
                Yz[kr] = c[0];
                Yz[ki] = c[1];
                x = Yz[kr];
                // Y[k] = S1[k] * H1[k] + j * S2[k] * H2[k]
                //      = reT[k] - imY[k] + j*(imT[k] + reY[k])
                Yz[kr] = tempz[kr] - Yz[ki];
                Yz[ki] = tempz[ki] + x;

                // handle conjugate-symmetric upper half of H1 and H2
                c = cmult(tempz[lr], tempz[li], H1z[kr], -H1z[ki]);
                tempz[lr] = c[0];
                tempz[li] = c[1];
                c = cmult(Yz[lr], Yz[li], H2z[kr], -H2z[ki]);
                Yz[lr] = c[0];
                Yz[li] = c[1];
                x = Yz[lr];
                // Y[k] = S1[k] * H1[k] + j * S2[k] * H2[k]
                //      = reT[k] - imY[k] + j*(imT[k] + reY[k])
                Yz[lr] = tempz[lr] - Yz[li];
                Yz[li] = tempz[li] + x;
            }
        } else {
            // DC and Nyquist frequency
            Yz[0] = tempz[0] * H1z[0];
            Yz[Nyr] = tempz[Nyr] * H1z[Nyr];

            // interleaved freq domain filter
            for (k = 1; k < Ny; k++) {
                int kr = 2 * k;
                int ki = kr + 1;
                l = Nz - k;
                int lr = 2 * l;
                int li = lr + 1;
                double[] c = cmult(tempz[kr], tempz[ki], H1z[kr], H1z[ki]);
                Yz[kr] = c[0];
                Yz[ki] = c[1];
                c = cmult(tempz[lr], tempz[li], H1z[kr], -H1z[ki]);
                Yz[lr] = c[0];
                Yz[li] = c[1];
            }
        }
        ifft(Yz, Wz, Nz);
        return 0L;
    }

    static void singleFir(double[] Wz, double[] X, double[] tempz, double[] Yz,
        double[] Hz, int N) {
        double[] Xz = new double[2 * N];
        for (int i = 0; i < X.length; i++) {
            Xz[2 * i] = X[i];
        }
        dualFir(Wz, Xz, tempz, Yz, Hz, Hz, N);
    }

    /**
     * Frequency domain filtering of one real sequence using half-size FFT.
     * 
     * @param Wz [N/4] twiddling factors (complex interleaved, of half-size FFT)
     * @param Xz [N] real input vector
     * @param tempz [N] temporary buffer
     * @param Y [N] real output vector
     * @param Hz [N/4 + 1] DFT of IR (complex interleaved)
     * @param N length of input and output vectors
     * @return
     */
    // TODO: test
    public static long fir(double[] Wz, double[] X, double[] tempz, double[] Y,
        double[] Hz, int N) {
        int Nz = N / 2;
        int Ny = Nz / 2;
        fftr(X, Wz, Nz);

        // DC and Nyquist frequency
        Y[0] = tempz[0] * Hz[0];
        Y[Ny] = tempz[Ny] * Hz[Ny];

        // interleaved freq domain filter
        for (int k = 1; k < Ny; k++) {
            int kr = 2 * k;
            int ki = kr + 1;
            int l = Nz - k;
            int lr = 2 * l;
            int li = lr + 1;
            double[] c = cmult(tempz[kr], tempz[ki], Hz[kr], Hz[ki]);
            Y[kr] = c[0];
            Y[ki] = c[1];
            c = cmult(tempz[lr], tempz[li], Hz[kr], -Hz[ki]);
            Y[lr] = c[0];
            Y[li] = c[1];
        }
        ifftr(Y, Wz, Nz);
        return 0L;
    }

    //////////// helper methods /////////////

    /**
     * print complex vector
     */
    static void cprintf(double[] az, int Nz) {
        int i;
        for (i = 0; i < Nz; i++) {
            System.out.print(String.format("%d %2.5f %s%2.5fi\n", i, az[2 * i],
                az[2 * i + 1] >= 0 ? "+" : "", az[2 * i + 1]));
        }
    }

    /**
     * print real vector
     * 
     * @param a
     * @param N
     */
    static void rprintf(double[] a, int N) {
        int i;
        for (i = 0; i < N; i++) {
            System.out.print(String.format("%d %2.5f\n", i, a[i]));
        }
    }

    /**
     * bit reversal algorithm
     * 
     * @param i
     * @param m
     * @return
     */
    static int bitrev(int i, int m) {
        int ir, k;
        ir = 0;
        int n = 1 << m;
        for (k = 1; k <= m; k++) {
            n >>= 1;
            ir = ir | ((i >> (k - 1)) & 1) << (m - k);
        }
        return ir;
    }

    /**
     * dual log
     * 
     * @param n
     * @return
     */
    static int ldint(int n) {
        int p = -1;
        int m = n;
        while (m != 0) {
            m = m >> 1;
            p++;
        }
        // exact power?
        if (1 << p == n) {
            return p;
        }
        return p + 1;
    }

    /**
     * permuation of signal bins
     * 
     * @param xz complex signal
     * @param nz complex element count
     */
    static void perm(double[] xz, int nz) {
        int i, ir;
        double temp;
        int m = ldint(nz);
        for (i = 0; i < nz; i++) {
            ir = bitrev(i, m);
            // only swap once
            if (ir > i) {
                // swap
                temp = xz[2 * ir];
                xz[2 * ir] = xz[2 * i];
                xz[2 * i] = temp;
                temp = xz[2 * ir + 1];
                xz[2 * ir + 1] = xz[2 * i + 1];
                xz[2 * i + 1] = temp;
            }
        }
    }

    /**
     * raw implementation of radix-2 decimation-in-time FFT
     * 
     * @param Xz complex input
     * @param Wz complex twiddle factors
     * @param idir direction +/-1
     * @param Nz complex FFT size
     */
    static void fft1d(double[] Xz, double[] Wz, int idir, int Nz) {
        int i, u, v, inca, incb, incn, j, k, ell;
        double[] z = {0, 0};

        incn = Nz;
        inca = 1;
        // in-place butterfly (Danielson-Lanczos)
        while (incn > 1) {
            incn /= 2;
            incb = 2 * inca;
            i = 0;
            for (j = 0; j < incn; j++) {
                k = 0;
                for (ell = 0; ell < inca; ell++) {
                    u = i + ell;
                    v = u + inca;
                    int kr = 2 * k;
                    int ki = kr + 1;
                    int vr = 2 * v;
                    int vi = vr + 1;
                    int ur = 2 * u;
                    int ui = ur + 1;
                    z[0] = Wz[kr] * Xz[vr] - Wz[ki] * idir * Xz[vi];
                    z[1] = Wz[kr] * Xz[vi] + Wz[ki] * idir * Xz[vr];
                    Xz[vr] = Xz[ur] - z[0];
                    Xz[vi] = Xz[ur + 1] - z[1];
                    Xz[ur] = Xz[ur] + z[0];
                    Xz[ui] = Xz[ui] + z[1];
                    k += incn;
                }
                i += incb;
            }
            inca = incb;
        }

        // multiply by 1/n for inverse transform
        if (idir < 0) {
            double norm = 1.0 / Nz;
            for (i = 0; i < 2 * Nz; i++) {
                Xz[i] = Xz[i] * norm;
            }
        }
    }

    //    /**
    //     * raw implementation of radix-2 decimation-in-time FFT for real sequences
    //     * 
    //     * @param Xz complex input / output (Note that the Nyquist frequency value
    //     *        is put into the imaginary part of DC)
    //     * @param Wz complex twiddle factors for FFT of size N/2
    //     * @param idir direction +/-1
    //     * @param N complex signal length (twice the FFT size)
    //     */
    //    // TODO: test
    //    static void fftr1d(double[] X, double[] Wz, int idir, int N) {
    //        double[] wz = {0, 0};
    //        int Nz = N / 2;
    //        double pi2n = 2 * Math.PI / N;
    //        double cospi2n = Math.cos(pi2n);
    //        double sinpi2n = Math.sin(pi2n);
    //        double[] u = {0, 0}, v = {0, 0};
    //        double[] Xz = X;
    //        // N/2 forward transform
    //        if (idir == 1) {
    //            perm(Xz, Nz);
    //            fft1d(Xz, Wz, idir, Nz);
    //        }
    //
    //        // DC and Nyquist freq
    //        u[0] = Xz[0];
    //        u[1] = Xz[1];
    //        Xz[0] = u[0] + u[1];
    //        Xz[1] = u[0] - u[1];
    //
    //        // nontrivial freqs
    //        for (int k = 1; k < N / 4; k++) {
    //            int kr = 2 * k;
    //            int ki = kr + 1;
    //            int lr = N / 2 - kr;
    //            int li = lr + 1;
    //            u = wz;
    //            wz[0] = cospi2n * u[0] - sinpi2n * u[1];
    //            wz[1] = cospi2n * u[1] + sinpi2n * u[0];
    //            u[0] = Xz[kr];
    //            u[1] = Xz[ki];
    //            v[0] = Xz[lr];
    //            v[1] = Xz[li];
    //            Xz[kr] = .5 * (u[0] + v[0] + wz[0] * idir * (u[1] + v[1]) + wz[1]
    //                * (u[0] - v[0]));
    //            Xz[ki] = .5 * (u[1] - v[1] - wz[0] * idir * (u[0] - v[0]) + wz[1]
    //                * (u[1] + v[1]));
    //            Xz[lr] = .5 * (u[0] + v[0] - wz[0] * idir * (u[1] + v[1]) - wz[1]
    //                * (u[0] - v[0]));
    //            Xz[li] = .5 * (-u[1] + v[1] - wz[0] * idir * (u[0] - v[0]) + wz[1]
    //                * (u[1] + v[1]));
    //        }
    //
    //        if (idir == -1) {
    //            Xz[0] *= .5;
    //            Xz[1] *= .5;
    //            // N/2 transform
    //            perm(Xz, Nz);
    //            fft1d(Xz, Wz, idir, Nz);
    //        }
    //
    //    }

    /**
     * forward real FFT with half length
     * 
     * @param X data [N] (real)
     * @param Yz output [N/2 + 1] (complex)
     * @param Wz [N/2] coefficients (complex interleaved, of full-size FFT)
     * @param N FFT size (real)
     */
    /*
    public static void fftr(double[] X, double[] Yz, double[] Wz, int N) {
        // complex interleaved sequence from real one
        int Nz = N / 2;
        double[] Xz = X;
        double[] X2z = new double[N];
        double[] W2z = new double[N / 2];

        // half-length twiddle factors
        for (int i = 0; i < W2z.length; i++) {
            W2z[i] = Wz[2 * i];
        }

        // permute input
        perm(Xz, Nz);

        // perform inverse transform
        fft1d(Xz, W2z, 1, Nz);

        int Nyz = Nz / 2;

        // symmetric FFT approach
        for (int k = 1; k < Nyz; k++) {
            int kr = 2 * k;
            int ki = kr + 1;
            int l = Nz - k;
            int lr = 2 * l;
            int li = lr + 1;
            // X2[k] = -j*(Z[k] - Z'[l]) / 2
            //      = (imZ[k] + imZ[l] + j*(reZ[l] - reZ[k])) / 2;
            X2z[kr] = (Xz[ki] + Xz[li]) * 0.5;
            X2z[ki] = (Xz[lr] - Xz[kr]) * 0.5;
            // X1[k] = (Z[k] + Z'[l]) / 2
            //      = (reZ[k] + reZ[l] + j*(imZ[k] - imZ[l])) / 2
            Xz[kr] = (Xz[kr] + Xz[lr]) * 0.5;
            Xz[ki] = (Xz[ki] - Xz[li]) * 0.5;

            // X[l] = reX[k] - j*imX[k];
            Xz[lr] = Xz[kr];
            Xz[li] = -Xz[ki];
            // X2[l] = reX2[k] - j*imX2[k];
            X2z[lr] = X2z[kr];
            X2z[li] = -X2z[ki];

            // Y[k] = X1[k] + W_2N^k * X2[k]
            // Y[2N - k] = X1[k] - W_2N^k * X2[k]
            double[] c = cmult(Wz[kr], Wz[ki], X2z[kr], X2z[ki]);
            Yz[kr] = Xz[kr] + c[0];
            c = cmult(Wz[kr], Wz[ki], X2z[kr], X2z[ki]);
            Yz[kr] = Xz[kr] + c[0];

        }
    }
    */

    static double[] cmult(double xr, double xi, double yr, double yi) {
        double[] z = {0, 0};
        double ac, bd;

        ac = xr * yr;
        bd = xi * yi;
        z[0] = ac - bd;
        z[1] = (xr + xi) * (yr + yi) - ac - bd;

        return z;
    }

    public static void main(String[] args) {
        int N = 32;
        int nz = N / 2;
        int i;
        double[] x1 = new double[N];
        double[] x2 = new double[N];
        double[] xz = new double[2 * N];
        double[] tempz = new double[2 * N];
        double[] wz = new double[N];
        //double[] hz = aa;
        double[] hz = Vectors.zeros(32);
        hz[0] = 1.;
        hz[2] = 1.;
        // create test signal
        for (i = 0; i < N; i++) {
            x1[i] = 0 - i;//Math.sin(Math.PI * i / 20) + Math.sin(Math.PI * i / 12);
            x2[i] = i;
        }
        
        // create FFT coefficients
        ffw(wz, N);

        // transform IR
        fft(hz, wz, nz);

        // interleave signals
        for (i = 0; i < N; i++) {
            xz[2 * i] = x1[i];
            xz[2 * i + 1] = x2[i];
        }
        System.out.println("*** original ***");
        FirFilter.cprintf(xz, N);

        dualFir(wz, xz, tempz, xz, hz, hz, nz);
        System.out.println("filtered xy\n");
        cprintf(xz, nz);
        System.out.println();
    }

    // example filter coeffs, fc = 0.05 lowpass
    static final double aa[] = {0.000952938986217, 0.001555347216713,
        0.002976368327833, 0.005135668181840, 0.008185514616811,
        0.012232925387171, 0.017313508506229, 0.023369751798990,
        0.030237910518926, 0.037646833816889, 0.045230583547699,
        0.052554733546986, 0.059154102941475, 0.064577737117313,
        0.068435556577005, 0.070440518911905, 0.070440518911905,
        0.068435556577005, 0.064577737117313, 0.059154102941475,
        0.052554733546986, 0.045230583547699, 0.037646833816889,
        0.030237910518926, 0.023369751798990, 0.017313508506229,
        0.012232925387171, 0.008185514616811, 0.005135668181840,
        0.002976368327833, 0.001555347216713, 0.000952938986217};
}
