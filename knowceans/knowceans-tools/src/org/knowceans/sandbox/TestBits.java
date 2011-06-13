/*
 * Created on Nov 10, 2003
 *
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
/*
 * Copyright (c) 2002 Gregor Heinrich.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESSED 
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY 
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/*
 * Created on Jan 9, 2005
 */
package org.knowceans.sandbox;


/**
 * TestBits sheds some light on Java unsigned number handling and conversion.
 * 
 * @author heinrich
 */
public class TestBits {
    public static void main(String[] args) {
        // start with 16 bit
        short x = -128;
        short x2 = -0x80;
        int rightY = 65408;
        
        int i1 = -784561419;
        long l1 = 3510405877l;              
                
        System.out.println(toBits(i1));
        System.out.println(toBits(l1));

        
        long l2 = i1 & 0x00000000ffffffffl;
        System.out.println(toBits(l2));
        
    }

    public static String toBits(short a) {
        String s = "";
        int p = 1;
        int x = 0;
        for (int i = 0; i < Short.SIZE; i++) {
            x = (a & p) / p;
            s = x + s;
            p *= 2;
        }

        return s;
    }

    public static String toBits(int a) {
        String s = "";
        int p = 1;
        int x = 0;
        for (int i = 0; i < Integer.SIZE; i++) {
            x = (a & p) / p;
            s = x + s;
            p *= 2;
        }

        return s;
    }

    public static String toBits(long a) {
        String s = "";
        long p = 1;
        long x = 0;
        for (int i = 0; i < Long.SIZE; i++) {
            x = (a & p) / p;
            s = x + s;
            p *= 2;
        }

        return s;
    }

    
    public static int toUnsigned(short a) {
        int b = a;
        if (a < 0)
            b = a + Short.MAX_VALUE + (Short.MIN_VALUE - a);
        return b;

    }
}
