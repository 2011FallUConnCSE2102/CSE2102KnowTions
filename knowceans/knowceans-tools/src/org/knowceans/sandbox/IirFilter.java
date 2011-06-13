/*
 * Created on Feb 28, 2010
 */
package org.knowceans.sandbox;

import java.util.Arrays;

/**
 * IirFilter implements a simple IIR filter using a cascade of biquad filters
 * (second-order sections).
 * 
 * @author gregor
 */
public class IirFilter {

    public IirFilter(int channels, int stages, int framesize) {
        v = new double[stages][channels * 3];
        init(channels, stages, framesize);
    }

    public void init(int channels, int stages, int framesize) {
        nchannel = channels;
        nstage = stages;
        nframe = framesize;
        for (int i = 0; i < nstage; i++) {
            Arrays.fill(v[i], 0);
        }
    }

    int nchannel;
    int nstage;
    int nframe;
    double[][] a;
    double[][] b;
    double[][] v;

    /**
     * perform processing for one frame
     * 
     * @param inframe
     * @param outframe
     */
    public void nextFrame(double[] inframe, double[] outframe, double gain) {
        // for each sample
        for (int i = 0; i < outframe.length; i++) {
            outframe[i] = inframe[i];
            // for each stage
            for (int k = 0; k < nstage; k++) {
                // update state using input and denominator coefficients
                v[k][0] = outframe[i] - a[k][0] * v[k][1] - a[k][1] * v[k][2];
                // update output using numerator coefficients and state
                outframe[i] = b[k][0] * v[k][0] + b[k][1] * v[k][1] + b[k][2]
                    * v[k][2];
                // progress in time
                v[k][2] = v[k][1];
                v[k][1] = v[k][0];
            }
            // global gain (may be factored into b[0])
            outframe[i] *= gain;
        }
    }

    public static void main(String[] args) {

        // a0 omitted: always 1
        double[][] aa = { {-1.867051864128537, 0.875313511144923},
            {-1.921144065558925, 0.946622554622648}};
        double[][] bb = { {1.0, 2.0, 1.0}, {1.0, 2.0, 1.0}};
        double gain = 1.241996358804925e-05;

        int N = 100;
        int i;
        double[] x = new double[N];
        double[] y = new double[N];
        // create test signal
        for (i = 0; i < N; i++) {
            x[i] = Math.sin(Math.PI * i / 20) + Math.sin(Math.PI * i / 12);
        }

        System.out.println("*** original ***");
        FirFilter.rprintf(x, N);

        // assemble filter
        IirFilter filter = new IirFilter(1, 2, 10);
        filter.a = aa;
        filter.b = bb;

        double[] inframe = new double[filter.nframe];
        double[] outframe = new double[filter.nframe];

        // process signal
        for (i = 0; i < N; i += filter.nframe) {
            System.arraycopy(x, i, inframe, 0, filter.nframe);
            filter.nextFrame(inframe, outframe, gain);
            System.arraycopy(outframe, 0, y, i, filter.nframe);
        }

        System.out.println("*** filtered ***");
        FirFilter.rprintf(y, N);
    }
}
