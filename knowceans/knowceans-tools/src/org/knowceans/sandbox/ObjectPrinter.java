/*
 * Created on Nov 15, 2003
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

package org.knowceans.sandbox;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * Class to analyse serialised objects. Analyses all fields, trying to access private
 * fields. Also tries to evaluate all methods with a non-void return parameter and
 * no arguments. The results of the field and method evaluation are displayed, and 
 * more complex classes displayed with their private and public fields. If the childclass
 * is a Collection descendent, the contained elements are displayed if non-null.
 * The class is rather a hack. Can be restructured to dump debug outputs at runtime.
 * 
 * @author heinrich
 */
public class ObjectPrinter {

    boolean overridePrivate = true;
    Object o = null;

    public static void main(String[] args) {
//		String z = "c:/incoming/ebiziao-workswithPublikationenSubdirUndPDFMask.test";
		String z = "c:/incoming/ebiziao-Pub-link!=DL.test";
		
        if (args.length > 0 && args[0] != null)
            z = args[0];

        ObjectPrinter s = new ObjectPrinter(new File(z));

		s.analyse();

    }

    public ObjectPrinter(File z) {

        try {
            ObjectInputStream in = new ObjectInputStream(new FileInputStream(z));
            o = in.readObject();
            in.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }

    public ObjectPrinter(Object c) {

        o = c;
    }

    public void analyse() {
        try {
            Object c = o;
            StringBuffer b = new StringBuffer();

            b.append(modifierString(c.getClass().getModifiers()) + c.getClass() + " {\n\n");
            b.append("// Fields\n");

            Field f[] = null;
			f = c.getClass().getFields();
			Set fx = new HashSet(Arrays.asList(f));
            if (overridePrivate)
                fx.addAll(Arrays.asList(c.getClass().getDeclaredFields()));
            f = (Field[]) fx.toArray(new Field[0]);
            

            for (int i = 0; i < f.length; i++) {
                f[i].setAccessible(overridePrivate);
                b.append(modifierString(f[i].getModifiers()) + f[i].getType().getName() + " " + f[i].getName() + "\n");

                Object x = f[i].get(c);
                b.append(printValue(x));
                if (x != null
                    && !f[i].getType().isPrimitive()
                    && !x.getClass().getName().startsWith("java.lang")
                    && !f[i].getType().isArray()) {

                    b.append(showSubObject(f[i].getName(), f[i].get(c)) + "\n");

                }
            }

            b.append("\n// Methods\n");

            Method m[] = c.getClass().getMethods();
            for (int i = 0; i < m.length; i++) {

                b.append(
                    modifierString(m[i].getModifiers()) + m[i].getReturnType().getName() + " " + m[i].getName() + "(");
                Class p[] = m[i].getParameterTypes();
                for (int j = 0; j < p.length; j++) {
                    b.append(p[j].getName());
                    if (j < p.length - 1)
                        b.append(", ");
                }
                b.append(");\n");
                if (!m[i].getReturnType().getName().equals("void") && p.length == 0) {
                    if (!m[i].getReturnType().isArray()
                        && !m[i].getReturnType().isPrimitive()
                        && !m[i].getReturnType().getName().startsWith("java.lang"))
                        b.append(showSubObject(m[i].getName() + "()", m[i].invoke(c, new Object[0])) + "\n");
                }
            }
            b.append("}");
            System.out.println(b);

        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        } catch (InvocationTargetException e) {
            e.printStackTrace();
        }
    }

    /**
     * @param object
     * @return
     */
    private String printValue(Object o) {
        StringBuffer b = new StringBuffer();
        if (o != null) {

            if (o.getClass().isArray()) {
                Class ca = o.getClass().getComponentType();
                if (ca.isPrimitive()) {

                    // Java <1.5 would make it necessary to element-wise
                    // copy the primitive array to an object type, thus direct printing.
                    int l = Array.getLength(o);
                    b.append("[");
                    for (int i = 0; i < Array.getLength(o); i++) {
                        b.append(Array.get(o, i));
                        if (i < l - 1)
                            b.append(", ");

                    }
                    b.append("];\n");
                } else {
                    // object arrays
                    b.append("  = " + Arrays.asList((Object[]) o) + ";\n");
                }
            } else if (Collection.class.isAssignableFrom(o.getClass())) {
                b.append("  = " + o + ";\n");
            } else {
                b.append("  = " + o + ";\n");
            }
        }
        return b.toString();
    }

    /**
     * @param i
     * @return
     */
    private String modifierString(int i) {
        StringBuffer a = new StringBuffer();
		if (Modifier.isPrivate(i))
			a.append("private ");
		if (Modifier.isProtected(i))
			a.append("protected ");
		if (Modifier.isPublic(i))
			a.append("public ");
		if (Modifier.isSynchronized(i))
			a.append("synchronized ");
        if (Modifier.isAbstract(i))
            a.append("abstract ");
        if (Modifier.isInterface(i))
            a.append("interface ");
        if (Modifier.isNative(i))
            a.append("native ");
        if (Modifier.isStatic(i))
            a.append("static ");
		if (Modifier.isFinal(i))
			a.append("final ");
        if (Modifier.isTransient(i))
            a.append("transient ");
        if (Modifier.isVolatile(i))
            a.append("volatile ");
        return a.toString();
    }

    public String showSubObject(String name, Object c) {
        if (c == null)
            return null;
        try {
            StringBuffer b = new StringBuffer();
            b.append("\n          _________________ " + name + " : _________________\n\n");
            b.append("          " + modifierString(c.getClass().getModifiers()) + c.getClass() + " {\n\n");

			Field f[] = null;
			f = c.getClass().getFields();
			Set fx = new HashSet(Arrays.asList(f));
			if (overridePrivate)
				fx.addAll(Arrays.asList(c.getClass().getDeclaredFields()));
			f = (Field[]) fx.toArray(new Field[0]);
            for (int i = 0; i < f.length; i++) {
                f[i].setAccessible(overridePrivate);
                b.append(
                    "          "
                        + modifierString(f[i].getModifiers())
                        + f[i].getType().getName()
                        + " "
                        + f[i].getName()
                        + "\n");
                Object x = f[i].get(c);
                b.append("          " + printValue(x));

                if (Collection.class.isAssignableFrom(c.getClass())
                    && f[i].getName().equals("elementData")
                    && !name.startsWith("ELEMENT ")) {
                    Object[] xa = (Object[]) x;
                    for (int j = 0; j < xa.length; j++) {
                        if (xa[j] != null)
                            b.append(showSubObject("ELEMENT " + j, xa[j]) + "\n");

                    }

                }
            }

            b.append("          }\n");
            b.append("          __________________________________________\n");
            return b.toString();

        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
            //        } catch (InvocationTargetException e) {
            //            e.printStackTrace();
        }
        return null;

    }

}
