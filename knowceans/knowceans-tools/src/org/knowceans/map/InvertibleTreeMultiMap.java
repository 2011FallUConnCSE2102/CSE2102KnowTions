/*
 * Copyright (c) 2005-2006 Gregor Heinrich. All rights reserved. Redistribution and
 * use in source and binary forms, with or without modification, are permitted
 * provided that the following conditions are met: 1. Redistributions of source
 * code must retain the above copyright notice, this list of conditions and the
 * following disclaimer. 2. Redistributions in binary form must reproduce the
 * above copyright notice, this list of conditions and the following disclaimer
 * in the documentation and/or other materials provided with the distribution.
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESSED OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package org.knowceans.map;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This object maps each key to a Set of values and holds each member of this
 * Set as a link to a Set of keys, i.e., depicts a IMultiMap (object -> set)
 * with an inverse IMultiMap. With getInverse(), the (forward) key set can
 * easily be found from the (forward) value. This effectively implements an M to
 * N relation.
 * <p>
 * In relational terms, this class implements an m:n relation.
 * <p>
 * 
 * @author heinrich
 */
public class InvertibleTreeMultiMap<X, Y> extends TreeMultiMap<X, Y>
// implements IInvertibleMultiMap<X, Y>
{

    /**
     * Comment for <code>serialVersionUID</code>
     */
    private static final long serialVersionUID = -4060547918794757464L;
    /**
     * keeps all inverse mappings, for each value that appears in the forward
     * values, one key is created that maps to the set of keys that point to it
     * in the forward map.
     */
    private HashMultiMap<Y, X> inverse;

    /**
     * some simple tests and demonstration for HashMultiMap
     * 
     * @param args
     */
    public static void main(String[] args) {
        InvertibleTreeMultiMap<String, Integer> m = new InvertibleTreeMultiMap<String, Integer>();
        m.add("c", 3);
        m.add("a", 1);
        m.add("b", 2);

        Set<Integer> d = new HashSet<Integer>();
        d.add(4);
        d.add(14);
        d.add(24);
        m.add("bb", 2);
        System.out.println(m);
        System.out.println(m.getInverse());

    }

    /**
     * 
     */
    public InvertibleTreeMultiMap() {
        super();
        inverse = new HashMultiMap<Y, X>();
    }

    /**
     * @param c
     */
    public InvertibleTreeMultiMap(Comparator< ? super X> c) {
        super(c);
        inverse = new HashMultiMap<Y, X>();
    }

    /**
     * @param m
     */
    public InvertibleTreeMultiMap(Map m) {
        super(m);
        inverse = new HashMultiMap<Y, X>();
    }

    /**
     * @param inverse allow to set an inverse map that could be a tree map with
     *        its own comparator.
     * @param m
     */
    public InvertibleTreeMultiMap(Map m, HashMultiMap<Y, X> inverse) {
        super(m);
        this.inverse = inverse;
    }

    @Override
    public void clear() {
        super.clear();
        inverse.clear();
    }

    /**
     * adds the value to the key's value set and the key to the value's key set.
     */
    public void add(X key, Y value) {
        super.add(key, value);
        inverse.add(value, key);
    }

    /**
     * put a new key-valueSet pair. Value must be a Set. The old set of values
     * must be removed if they don't , which places in the values set as well
     * for every entry of the value, a In the inverse map, the (forward) key
     * needs to be removed from the Set that is pointed to by its old value, and
     * the key put to the Set for the new value.
     */
    public Set<Y> put(X key, Set<Y> valSet) {
        // Put into the forward map, which is super.
        // returns a Set of old relations to the key,
        // which must be cleaned from the reverse side
        Set<Y> oldVal = super.put(key, valSet);

        // Update the reverse map. oldVal no longer maps to key.
        // if nothing was set for the key oldVal=null
        if (oldVal != null) {
            for (Y el : oldVal) {
                inverse.remove(el, key);
                if (inverse.get(el).size() == 0)
                    inverse.remove(el);
            }
        }
        // put in the new reverse relation of every value from
        // the set to the
        for (Y el : valSet) {
            inverse.add(el, key);
        }

        // Return the old value.
        return oldVal;
    }

    /**
     * removes the (forward) key and its elements from it value set from the map
     * that map to key.
     */
    @SuppressWarnings("unchecked")
    public Set<Y> remove(Object key) {
        // when subclassing, this results in a name clash because Map.remove
        // uses Object as key type.
        // public Set<Y> remove(K key) {
        X xkey = (X) key;
        // Remove the inverse mapping and return the value mapped by key.
        Set<Y> val = super.remove(key);
        // Update the reverse map. oldVal no longer maps to key.
        for (Y el : val) {
            inverse.remove(el, xkey);
            Set<X> inv = inverse.get(el);
            if (inv != null && inv.size() == 0)
                inverse.remove(el);
        }

        return val;
    }

    /**
     * removes the value from the key's mapping.
     */
    public void remove(X key, Y value) {
        super.remove(key, value);
        inverse.remove(value, key);
        if (inverse.get(value).size() == 0)
            inverse.remove(value);

    }

    /**
     * gets keys for a value as a Set.
     * 
     * @param val
     * @return
     */
    public Set<X> getInverseValue(Y val) {
        return inverse.get(val);
    }

    /**
     * returns the keys that match a given set of value elements
     * 
     * @param value
     * @return
     */
    public Set<X> getInverse(Set<Y> value) {
        return inverse.get(value);
    }

    /**
     * returns the hash map of the inverse (not that this is not a tree map).
     * 
     * @return
     */
    public IMultiMap<Y, X> getInverse() {
        return inverse;
    }

}
