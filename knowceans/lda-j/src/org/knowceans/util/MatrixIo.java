/*
 * Created on Aug 1, 2005
 */
package org.knowceans.util;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

/**
 * MatrixIo loads and saves a binary matrix in a zipped file
 *
 * @author heinrich
 */
public class MatrixIo {
    
    public static void main(String[] args) {
        
        double[][] a = {{2, 4},{1, 3},{7, 4}};
        
        saveBinaryMatrix("test.zip", a);
        double[][] b = loadBinaryMatrix("test.zip");
        
        for (int i = 0; i < b.length; i++) {
            for (int j = 0; j < b[i].length; j++) {
                System.out.print(" " + b[i][j]);
            }
            System.out.println();
        }
    }
    
    /**
     * loads a matrix from a binary file
     * 
     * @param filename
     * @return
     */
    public static double[][] loadBinaryMatrix(String filename) {
        int m, n;
        double[][] a = null;
        int i = 0, j = 0;
        try {

            DataInputStream dis = null;

            if (filename.endsWith(".zip")) {

                ZipFile f = new ZipFile(filename);
                String name = new File(filename).getName();
                dis = new DataInputStream(new BufferedInputStream(f
                    .getInputStream(f.getEntry(name.substring(0,
                        name.length() - 3)
                        + "bin"))));
            } else {
                dis = new DataInputStream(new BufferedInputStream(
                    new FileInputStream(filename)));
            }
            m = dis.readInt();
            n = dis.readInt();
            a = new double[m][n];
            for (i = 0; i < m; i++) {
                for (j = 0; j < n; j++) {
                    a[i][j] = dis.readFloat();
                }
            }
            dis.close();

        } catch (IOException e) {
            System.err.println(i + " " + j);
            e.printStackTrace();
        }
        return a;
    }

    /**
     * writes matrix to binary file. If the file name ends with zip
     * 
     * @param filename
     * @param a
     */
    public static void saveBinaryMatrix(String filename, double[][] a) {
        int i = 0, j = 0;

        DataOutputStream dis = null;
        try {
            if (filename.endsWith(".zip")) {
                ZipOutputStream zip = new ZipOutputStream(new FileOutputStream(
                    filename));
                String name = new File(filename).getName();
                zip.putNextEntry(new ZipEntry(name.substring(0,
                    name.length() - 3)
                    + "bin"));
                dis = new DataOutputStream(new BufferedOutputStream(zip));
            } else {
                dis = new DataOutputStream(new BufferedOutputStream(
                    new FileOutputStream(filename)));
            }
            dis.writeInt(a.length);
            dis.writeInt(a[0].length);
            for (i = 0; i < a.length; i++) {
                for (j = 0; j < a[0].length; j++) {
                    dis.writeFloat((float) a[i][j]);
                }
            }
            dis.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            System.err.println(i + " " + j);
            e.printStackTrace();
        }
    }
}
