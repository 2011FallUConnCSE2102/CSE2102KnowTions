/*
 * Created on 05.05.2006
 */
package org.knowceans.sandbox;

public class IntByteArray {

    public static void main(String[] args) {

        int u = -240;
        
        
        byte[] a = intToByteArray(u);
        int b = byteArrayToInt(a);

        System.out.println(u + " " + b);

    }

    public static byte[] intToByteArray(int value) {
        byte[] b = new byte[4];
        for (int i = 0; i < 4; i++) {
            int offset = (b.length - 1 - i) * 8;
            b[i] = (byte) ((value >>> offset) & 0xFF);
        }
        return b;
    }

    /**
     * Convert the byte array to an int.
     * 
     * @param b The byte array
     * @return The integer
     */
    public static int byteArrayToInt(byte[] b) {
        return byteArrayToInt(b, 0);
    }

    /**
     * Convert the byte array to an int starting from the given offset.
     * 
     * @param b The byte array
     * @param offset The array offset
     * @return The integer
     */
    public static int byteArrayToInt(byte[] b, int offset) {
        int value = 0;
        for (int i = 0; i < 4; i++) {
            int shift = (4 - 1 - i) * 8;
            value += (b[i + offset] & 0x000000FF) << shift;
        }
        return value;
    }

    static void printBits(int a) {
        System.out.println(Integer.toBinaryString(a));
    }

    static void printBits(byte[] b) {
        for (int i = 0; i < b.length; i++) {
            int u = 0xff;
            for (int j = 0; j < 8; j++) {
                System.out.print(b[i] & u);
            }
        }
    }
}
