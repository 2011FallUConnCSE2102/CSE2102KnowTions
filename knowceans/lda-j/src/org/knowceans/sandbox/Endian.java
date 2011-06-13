/*
 * Created on Apr 4, 2005
 */
package org.knowceans.sandbox;

/**
 * Endian converts byte values from big to little endian and vice versa.
 * 
 * @author heinrich
 */
public class Endian {

    public static final short shortL2B(short w) {
        return (short) ((w & 0x00ff) << 8 | (w & 0xff00));
    }

    public static final char charL2B(char w) {
        return (char) ((w & 0x00ff) << 8 | (w & 0xff00));
    }

    public static final int intL2B(int w) {
        return (w & 0x000000ff) << 24 | (w & 0x0000ff00) << 8
            | (w & 0x00ff0000) >> 8 | (w & 0xff000000) >> 24;
    }

    public static final long longL2B(long w) {
        return (w & 0x000000ff) << 24 | (w & 0x0000ff00) << 8
            | (w & 0x00ff0000) >> 8 | (w & 0xff000000) >> 24;
    }
    
    public static void main(String[] args) {
        //int a = 0x0330ffff;
        int a = 1025;
        int b = intL2B(a);
        System.out.println(a + " -> " + b);
        System.out.println("0x" + Integer.toHexString(a) + " -> 0x" + Integer.toHexString(b));

    }
}