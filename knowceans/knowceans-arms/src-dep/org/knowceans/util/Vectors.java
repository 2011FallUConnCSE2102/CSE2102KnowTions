/*
 * Created on Jun 18, 2004 To change the template for this generated file go to
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

import java.util.Vector;

/**
 * Static vector manipulation routines for Matlab porting. The routines work for
 * int and double partly; the class is extended as needed. The question remains
 * whether it makes sense to have a formal IntVector, analoguous to IntMatrix
 * that allows performing search and indexing operations, such as views, etc.
 * 
 * @author heinrich
 */
public class Vectors {

    /**
     * @param start
     * @param end
     * @param step
     * @return [start : step : end]
     */
    public static int[] range(int start, int end, int step) {

        int[] out = new int[(int) Math.floor((end - start) / step) + 1];
        for (int i = 0; i < out.length; i++) {
            out[i] = start + step * i;
        }
        return out;
    }

    /**
     * @param start
     * @param end
     * @return [start : end]
     */
    public static int[] range(int start, int end) {
        return range(start, end, end - start > 0 ? 1 : -1);
    }

    /**
     * create sequence [start : step : end] of double values. TODO: check
     * precision.
     * 
     * @param start
     *            double value of start, if integer, use "1.0" notation.
     * @param end
     *            double value of end, if integer, use "1.0" notation.
     * @param step
     *            double value of step size
     * @return
     */
    public static double[] range(double start, double end, double step) {

        double[] out = new double[(int) Math.floor((end - start) / step) + 1];
        for (int i = 0; i < out.length; i++) {
            out[i] = start + step * i;
            System.out.println(step * i + " " + i);
        }
        return out;
    }

    /**
     * @param start
     * @param end
     * @return [start : end]
     */
    public static double[] range(double start, double end) {
        return range(start, end, end - start > 0 ? 1 : -1);
    }

    /**
     * sum the elements of vec
     * 
     * @param vec
     * @return
     */
    public static double sum(double[] vec) {
        double sum = 0;
        for (int i = 0; i < vec.length; i++) {
            sum += vec[i];
        }
        return sum;

    }

    /**
     * cumulative sum of the elements, starting at element 0.
     * 
     * @param vec
     * @return vector containing the cumulative sum of the elements of vec
     */
    public static double[] cumsum(double[] vec) {
        double[] x = new double[vec.length];
        x[0] = vec[0];
        for (int i = 1; i < vec.length; i++) {
            x[i] = vec[i] + x[i - 1];
        }
        return x;
    }

    /**
     * maximum value in vec
     * 
     * @param vec
     * @return
     */
    public static int max(int[] vec) {
        int max = vec[0];
        for (int i = 1; i < vec.length; i++) {
            if (vec[i] > max)
                max = vec[i];
        }
        return max;
    }

    /**
     * maximum value in vec
     * 
     * @param vec
     * @return
     */
    public static double max(double[] vec) {
        double max = vec[0];
        for (int i = 1; i < vec.length; i++) {
            if (vec[i] > max)
                max = vec[i];
        }
        return max;
    }

    /**
     * minimum value in vec
     * 
     * @param vec
     * @return
     */
    public static int min(int[] vec) {
        int min = vec[0];
        for (int i = 1; i < vec.length; i++) {
            if (vec[i] < min)
                min = vec[i];
        }
        return min;
    }

    /**
     * minimum value in vec
     * 
     * @param vec
     * @return
     */
    public static double min(double[] vec) {
        double min = vec[0];
        for (int i = 1; i < vec.length; i++) {
            if (vec[i] < min)
                min = vec[i];
        }
        return min;
    }

    /**
     * @param x
     * @param y
     * @return [x y]
     */
    public static double[] concat(double[] x, double[] y) {
        double[] z = new double[x.length + y.length];
        System.arraycopy(x, 0, z, 0, x.length);
        System.arraycopy(y, 0, z, x.length, y.length);
        return z;
    }

    /**
     * w = [x y z]
     * 
     * @param x
     * @param y
     * @return [x y z]
     */
    public static double[] concat(double[] x, double[] y, double[] z) {
        double[] w = new double[x.length + y.length + z.length];
        System.arraycopy(x, 0, w, 0, x.length);
        System.arraycopy(y, 0, w, x.length, y.length);
        System.arraycopy(y, 0, w, x.length + y.length, z.length);
        return w;
    }

    /**
     * prints a double representation of the vector.
     * 
     * @param x
     * @return
     */
    public static String print(double[] x) {
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < x.length - 1; i++) {
            b.append(x[i]).append(" ");
        }
        b.append(x[x.length - 1]);
        return b.toString();
    }

    /**
     * prints a double representation of an array.
     * 
     * @param x
     * @return
     */
    public static String print(double[][] x) {
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < x.length - 1; i++) {
            b.append(print(x[i])).append("\n");
        }
        b.append(print(x[x.length - 1]));
        return b.toString();
    }

    /**
     * prints a double representation of the vector.
     * 
     * @param x
     * @return
     */
    public static String print(int[] x) {
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < x.length - 1; i++) {
            b.append(x[i]).append(" ");
        }
        b.append(x[x.length - 1]);
        return b.toString();
    }

    /**
     * prints a double representation of an array.
     * 
     * @param x
     * @return
     */
    public static String print(int[][] x) {
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < x.length - 1; i++) {
            b.append(print(x[i])).append("\n");
        }
        b.append(print(x[x.length - 1]));
        return b.toString();
    }

    /**
     * @param len
     * @param factor
     * @return factor * ones(1, len);
     */
    public static double[] ones(int len, double factor) {
        double[] x = new double[len];
        for (int i = 0; i < x.length; i++) {
            x[i] = 1;
        }
        return x;
    }

    /**
     * @param len
     * @param factor
     * @return factor * ones(1, len);
     */
    public static int[] ones(int len, int factor) {
        int[] x = new int[len];
        for (int i = 0; i < x.length; i++) {
            x[i] = factor;
        }
        return x;
    }

    /**
     * @param len
     * @return zeros(1, len)
     */
    public static double[] zeros(int len) {
        return new double[len];
    }

    /**
     * @param len
     * @return ones(1, len)
     */
    public static int[] ones(int len) {
        return ones(len, 1);
    }

    /**
     * cast a double[] to an int[]
     * 
     * @param vec
     * @return
     */
    public static int[] cast(double[] vec) {
        int[] ivec = new int[vec.length];
        for (int i = 0; i < ivec.length; i++) {
            ivec[i] = (int) vec[i];
        }
        return ivec;
    }

    /**
     * cast a double[] to an int[]
     * 
     * @param vec
     * @return
     */
    public static double[] cast(int[] vec) {
        double[] dvec = new double[vec.length];
        for (int i = 0; i < dvec.length; i++) {
            dvec[i] = (double) vec[i];
        }
        return dvec;
    }

    /**
     * find indices with val
     * 
     * @param vec
     * @param val
     * @return vector with 0-based indices.
     */
    public static int[] find(int[] vec, int val) {
        Vector v = new Vector();
        for (int i = 0; i < vec.length; i++) {
            if (vec[i] == val) {
                v.add(new Integer(i));
            }
        }
        int[] vv = new int[v.size()];
        for (int i = 0; i < vv.length; i++) {
            vv[i] = ((Integer) v.get(i)).intValue();
        }
        return vv;
    }

    /**
     * returns a copy of the vector elements with the given indices in the
     * original vector.
     * 
     * @param indices
     * @return
     */
    public static double[] subVector(double[] vec, int[] indices) {
        double[] x = new double[indices.length];
        for (int i = 0; i < x.length; i++) {
            x[i] = vec[indices[i]];
        }
        return x;
    }

    /**
     * returns a copy of the vector elements with the given indices in the
     * original vector.
     * 
     * @param cols
     * @return
     */
    public static int[] subVector(int[] vec, int[] indices) {
        int[] x = new int[indices.length];
        for (int i = 0; i < x.length; i++) {
            x[i] = vec[indices[i]];
        }
        return x;
    }

    /**
     * @param weights
     * @param i
     * @param j
     * @return
     */
    public static double[] subVector(double[] vec, int start, int end) {
        double[] x = new double[end - start + 1];
        for (int i = 0; i <= end - start; i++) {
            x[i] = vec[start + i];
        }
        return x;
    }

    /**
     * set the elements of vec at indices with the respective replacements.
     * TODO: implement views as in the colt library
     * 
     * @param vec
     * @param indices
     * @param replacements
     * @return
     */
    public static void setSubVector(int[] vec, int[] indices, int[] replacements) {
        for (int i = 0; i < indices.length; i++) {
            vec[indices[i]] = replacements[i];
        }
    }

    /**
     * set the elements of vec at indices with the replacement. TODO: implement
     * views as in the colt library
     * 
     * @param vec
     * @param indices
     * @param replacement
     * @return
     */
    public static void setSubVector(int[] vec, int[] indices, int replacement) {
        for (int i = 0; i < indices.length; i++) {
            vec[indices[i]] = replacement;
        }
    }

    /**
     * add a scalar to the vector
     * 
     * @param vec
     * @param scalar
     */
    public static void add(int[] vec, int scalar) {
        for (int i = 0; i < vec.length; i++) {
            vec[i] += scalar;
        }
    }

    /**
     * set the elements of a copy of vec at indices with the respective
     * replacements. TODO: implement views as in the colt library
     * 
     * @param vec
     * @param indices
     * @param replacements
     * @return the copied vector with the replacements;
     */
    public static int[] setSubVectorCopy(int[] vec, int[] indices,
        int[] replacements) {
        int[] x = new int[vec.length];
        for (int i = 0; i < indices.length; i++) {
            x[indices[i]] = replacements[i];
        }
        return x;
    }

    /**
     * copies a the source to the destination
     * 
     * @param alpha
     * @return
     */
    public static double[] copy(double[] source) {
        if (source == null)
            return null;
        double[] dest = new double[source.length];
        System.arraycopy(source, 0, dest, 0, source.length);
        return dest;
    }

    /**
     * copies a the source to the destination
     * 
     * @param alpha
     * @return
     */
    public static int[] copy(int[] source) {
        if (source == null)
            return null;
        int[] dest = new int[source.length];
        System.arraycopy(source, 0, dest, 0, source.length);
        return dest;
    }

    /**
     * multiplicates the vector with a scalar. The argument is modified.
     * 
     * @param ds
     * @param d
     * @return
     */
    public static double[] mult(double[] ds, double d) {
        for (int i = 0; i < ds.length; i++) {
            ds[i] *= d;
        }
        return ds;

    }

    /**
     * multiplicates the vector with a vector (inner product). The argument is
     * not modified.
     * 
     * @param ds
     * @param d
     * @return
     */
    public static double mult(double[] ds, double[] dt) {
        if (ds.length != dt.length)
            throw new IllegalArgumentException("Vector dimensions must agree.");
        double s = 0;
        for (int i = 0; i < ds.length; i++) {
            s += ds[i] * dt[i];
        }
        return s;

    }    
}
