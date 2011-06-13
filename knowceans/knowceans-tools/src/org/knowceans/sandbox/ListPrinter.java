/*
 * Created on Nov 11, 2003
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

import java.util.Collection;
import java.util.Iterator;
import java.util.Vector;

/**
 * Convienience class to print lists. Should also work with
 * other collections. Users can set the public fields entry* and list*
 * to control printing behaviour. Further, a simple line breaking
 * mechanism is implemented.
 * 
 * @author heinrich
 */
public class ListPrinter<E> extends Vector<E> {

    /**
     * Comment for <code>serialVersionUID</code>
     */
    private static final long serialVersionUID = 3256438097242765105L;
    private int linewidth = 0;
    private int count = 0;
    public String entryStart = "'";
    public String entryEnd = "'";
    public String entrySeparator = ", ";
    public String listStart = "[";
    public String listEnd = "]";

    /**
     * empty printer (backed by a Vector).
     */
    public ListPrinter() {
        super();
    }
    /**
     * printer for Collection.
     * @param c
     */
    public ListPrinter(Collection<E> c) {
        super(c);
    }

    /**
     * empty printer (backed by a Vector) with
     * an approximate line width.
     * @param linewidth
     */
    public ListPrinter(int linewidth) {
        super();
        this.linewidth = linewidth;
    }

    /**
     * create an empty ListPrinter object
     * from the configuration of an existing one.
     * @param lp
     */
    public ListPrinter(ListPrinter< ? extends E> lp) {
        super();
        this.entryStart = lp.entryStart;
        this.entrySeparator = lp.entrySeparator;
        this.entryEnd = lp.entryEnd;
        this.linewidth = lp.linewidth;
        this.listStart = lp.listStart;
        this.listEnd = lp.listEnd;
    }

    /**
     * create ListPrinter object for the Collection
     * from the configuration of an existing ListPrinter.
     * @param lp
     */
    public ListPrinter(ListPrinter lp, Collection<E> c) {
        super(c);
        this.entryStart = lp.entryStart;
        this.entrySeparator = lp.entrySeparator;
        this.entryEnd = lp.entryEnd;
        this.linewidth = lp.linewidth;
        this.listStart = lp.listStart;
        this.listEnd = lp.listEnd;
    }

    /**
     * resets the printer content, but keeps the 
     * configuration. 
     * @param c
     */
    public void set(Collection<E> c) {
        super.clear();
        super.addAll(c);
    }
    /**
     * printer for Collection with an approximate
     * line width.
     * @param linewidth
     * @param c
     */
    public ListPrinter(int linewidth, Collection<E> c) {
        super(c);
        this.linewidth = linewidth;
    }

    /**
     * worker method to print out the list.
     */
    public synchronized String toString() {
        StringBuffer buf = new StringBuffer();
        buf.append(listStart);

        Iterator i = iterator();
        boolean hasNext = i.hasNext();
        while (hasNext) {
            Object o = i.next();
            String x = (o == this) ? "(this Collection)" : String.valueOf(o);
            count += x.length()
                + entryStart.length()
                + entryEnd.length()
                + entrySeparator.length();
            if (linewidth > 0 && count > linewidth) {
                buf.append("\n");
                count = 0;
            }
            buf.append(entryStart);
            buf.append(x);
            buf.append(entryEnd);
            hasNext = i.hasNext();
            if (hasNext)
                buf.append(entrySeparator);

        }
        buf.append(listEnd);
        return buf.toString();
    }
}