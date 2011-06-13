package org.knowceans.sandbox;

import java.util.Arrays;

/**
 * Convolver performs fast convolution on the input signals, using the FIR
 * filter on a signal vector.
 * 
 * @author gregor
 */
// TODO: test
public class Convolver {

    double[] save, HHz, Wz, tempz, XX;

    /**
     * performs overlap save convolution on the input signals xx
     * 
     * @param xx input signal frame
     * @param yy output signal frame (length = input + IR - 1)
     * @param hh time-domain impulse response
     * @param framesize 2^n
     */
    public static void convolve(double[] xx, double[] yy, double[] hh) {
        int i, j;
        Arrays.fill(yy, 0);
        // for all signal samples
        for (i = 0; i < xx.length; i++) {
            // for all impulse response samples
            for (j = 0; j < hh.length; j++) {
                yy[i + j] += xx[i] * hh[j];
            }
        }
    }

    /**
     * @param hh
     * @param framesize
     */
    public void init(double[] hh) {

        // double the size of IR length next power of 2
        int size = 2 << FirFilter.ldint(hh.length);
        // double fft size
        HHz = new double[2 * size];
        for (int i = 0; i < hh.length; i++) {
            HHz[2 * i] = hh[i];
        }
        Wz = new double[size];
        FirFilter.ffw(Wz, size);
        FirFilter.fft(HHz, Wz, size);
        tempz = new double[2 * size];
        save = new double[size / 2];
        XX = new double[size];

    }

    /**
     * performs overlap save convolution on the input signals xx
     * 
     * @param xx input signal frame vector
     * @param yy output signal frame vector
     * @param hh time-domain impulse response
     * @param framesize 2^n
     */
    public void overlapSave(double[] xx, double[] yy) {
        // prepend saved frame
        System.arraycopy(save, 0, XX, 0, save.length);
        System.arraycopy(xx, 0, XX, save.length, xx.length);
        System.arraycopy(xx, 0, save, 0, xx.length);
        // filter one frame
        FirFilter.singleFir(Wz, XX, tempz, XX, HHz, xx.length);
        // discard first half of output
        System.arraycopy(XX, xx.length, yy, 0, xx.length);
    }

}
