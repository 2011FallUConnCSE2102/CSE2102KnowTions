/*
 * Created on Apr 4, 2005
 */
package org.knowceans.sandbox;

import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;

/**
 * TestBinOutput
 * 
 * @author heinrich
 */
public class TestBinOutput {

    public static void main(String[] args) {
        String filename = "test.bin";
        try {
//            BufferedOutputStream bos = new BufferedOutputStream(
//                new FileOutputStream(filename));
            DataOutputStream bw = new DataOutputStream(new FileOutputStream(filename));
            double x = Math.PI;
            int y = 4;
            bw.writeDouble(x);
            bw.writeInt(y);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        double xx = 0;
        int yy = 0;
        
        try {
            DataInputStream dis = new DataInputStream(new FileInputStream(filename));
            xx = dis.readDouble();
            yy = dis.readInt();
            
        } catch (FileNotFoundException e1) {
            e1.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(xx + " " + yy);
        

    }

}
