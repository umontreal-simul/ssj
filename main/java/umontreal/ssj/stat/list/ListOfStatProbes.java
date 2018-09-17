/*
 * Class:        ListOfStatProbes
 * Description:  List of statistical probes
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist 
 * @since        2007
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.stat.list;

import umontreal.ssj.util.PrintfFormat;
import java.util.List;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.RandomAccess;
import umontreal.ssj.stat.StatProbe;

/**
 * Represents a list of statistical probes that can be managed
 * simultaneously. Each element of this list is a
 * @ref umontreal.ssj.stat.StatProbe instance which can be obtained and
 * manipulated.
 *
 * When constructing a list of statistical probes, one specifies the concrete
 * subclass of the  @ref umontreal.ssj.stat.StatProbe objects in it. One then
 * creates an empty list of probes, and fills it with statistical probes. If
 * the list is not intended to be modified, one can then use the
 * #setUnmodifiable to prevent any change in the contents of the list.
 *
 * Each list of statistical probes can have a global name describing the
 * contents of its elements, and local names associated with each individual
 * probe. For example, a list of statistical probes for the waiting times can
 * have the global name `Waiting times` while the individual probes have
 * local names `type 1`, `type 2`, etc. These names are used for formatting
 * reports.
 *
 * Facilities are provided to fill arrays with sums, averages, etc. obtained
 * from the individual statistical probes. Methods are also provided to
 * manipulate the contents of the list. However, one should always call
 * `init` immediately after adding or removing statistical probes in the
 * list.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ListOfStatProbes<E extends StatProbe>
                                implements Cloneable, List<E>, RandomAccess {

   // probes must implement RandomAccess, otherwise this class must not implement RandomAccess.
   private List<E> probes;
   private List<ArrayOfObservationListener> listeners = new ArrayList<ArrayOfObservationListener>();
   protected boolean collect = true;
   protected boolean broadcast = false;
   protected String name;

   /**
    * Constructs an empty list of statistical probes.
    */
   public ListOfStatProbes() {
      probes = new ArrayList<E>();
   }

   /**
    * Constructs an empty list of statistical probes with name `name`.
    *  @param name         the name of the new list.
    */
   public ListOfStatProbes (String name) {
      probes = new ArrayList<E>();
      this.name = name;
   }

   /**
    * Returns the global name of this list of statistical probes.
    *  @return the global name of the list.
    */
   public String getName() {
      return name;
   }

   /**
    * Sets the global name of this list to `name`.
    *  @param name         the new global name of the list.
    */
   public void setName (String name) {
      this.name = name;
   }

   /**
    * Determines if this list of statistical probes is modifiable, i.e.,
    * if probes can be added or removed. Any list of statistical probes is
    * modifiable by default, until one calls the  #setUnmodifiable method.
    */
   public boolean isModifiable() {
      return probes instanceof ArrayList;
   }

   /**
    * Forbids any future modification to this list of statistical probes.
    * After this method is called, any attempt to modify the list results
    * in an exception. Setting a list unmodifiable can be useful if some
    * data structures are defined depending on the probes in the list.
    */
   public void setUnmodifiable() {
      if (isModifiable())
         probes = Collections.unmodifiableList (probes);
   }

   /**
    * Initializes this list of statistical probes by calling
    * umontreal.ssj.stat.StatProbe.init on each element.
    */
   public void init() {
      for (StatProbe probe : probes) {
         if (probe != null)
            probe.init();
      }
   }

   /**
    * For each probe in the list, computes the sum by calling
    * umontreal.ssj.stat.StatProbe.sum, and stores the results into the
    * array `s`. This method throws an exception if the size of `s`
    * mismatches with the size of the list.
    *  @param s            the array to be filled with sums.
    *  @exception NullPointerException if `s` is `null`.
    *  @exception IllegalArgumentException if `s.length` does not
    * correspond to  #size.
    */
   public void sum (double[] s) {
      if (s.length != size())
         throw new IllegalArgumentException
            ("Invalid length of the given array: given length is " + s.length +
             ", required length is " + size());
      int i = 0;
      for (StatProbe probe : probes)
         s[i++] = probe == null ? Double.NaN : probe.sum();
   }

   /**
    * For each probe in this list, computes the average by calling
    * umontreal.ssj.stat.StatProbe.average, and stores the results into
    * the array `a`. This method throws an exception if the size of `s`
    * mismatches with the size of the list.
    *  @param a            the array to be filled with averages.
    *  @exception NullPointerException if `a` is `null`.
    *  @exception IllegalArgumentException if `a.length` does not
    * correspond to  #size.
    */
   public void average (double[] a) {
      if (a.length != size())
         throw new IllegalArgumentException
            ("Invalid length of the given array: given length is " + a.length +
             ", required length is " + size());
      int i = 0;
      for (StatProbe probe : probes)
         a[i++] = probe == null ? Double.NaN : probe.average();
   }

   /**
    * Determines if this list of statistical probes is collecting values.
    * Each probe of the list could or could not be collecting values. The
    * default is `true`.
    *  @return the status of statistical collecting.
    */
   public boolean isCollecting() {
      return collect;
   }

   /**
    * Sets the status of the statistical collecting mechanism to `c`. A
    * `true` value turns statistical collecting ON, a `false` value turns
    * it OFF.
    *  @param c            the status of statistical collecting.
    */
   public void setCollecting (boolean c) {
      collect = c;
   }

   /**
    * Determines if this list of statistical probes is broadcasting
    * observations to registered observers. The default is `false`.
    *  @return the status of broadcasting.
    */
   public boolean isBroadcasting() {
      return broadcast;
   }

   /**
    * Sets the status of the observation broadcasting mechanism to `b`. A
    * `true` value turns broadcasting ON, a `false` value turns it OFF.
    *  @param b            the status of broadcasting.
    */
   public void setBroadcasting (boolean b) {
      broadcast = b;
   }

   /**
    * Adds the observation listener `l` to the list of observers of this
    * list of statistical probes.
    *  @param l            the new observation listener.
    *  @exception NullPointerException if `l` is `null`.
    */
   public void addArrayOfObservationListener (ArrayOfObservationListener l) {
      if (l == null)
         throw new NullPointerException();
      if (!listeners.contains (l))
         listeners.add (l);
   }

   /**
    * Removes the observation listener `l` from the list of observers of
    * this list of statistical probes.
    *  @param l            the observation listener to be deleted.
    */
   public void removeArrayOfObservationListener (ArrayOfObservationListener l) {
      listeners.remove (l);
   }

   /**
    * Removes all observation listeners from the list of observers of this
    * list of statistical probes.
    */
   public void clearArrayOfObservationListeners() {
      listeners.clear();
   }

   /**
    * Notifies the observation `x` to all registered observers if
    * broadcasting is ON. Otherwise, does nothing.
    */
   public void notifyListeners (double[] x) {
      if (!broadcast)
         return;
      // We could also use the enhanced for loop here, but this is less efficient.
      final int nl = listeners.size();
      for (int i = 0; i < nl; i++)
         listeners.get (i).newArrayOfObservations (this, x);
   }

   /**
    * Formats a report for each probe in the list of statistical probes.
    * The returned string is constructed by using `StatProbe.report
    * (getName(), this)`.
    *  @return the report formatted as a string.
    */
   public String report() {
      return StatProbe.report (name, this);
   }

   /**
    * Clones this object. This makes a shallow copy of this list, i.e.,
    * this does not clone all the probes in the list. The created clone is
    * modifiable, even if the original list is unmodifiable.
    */
   public ListOfStatProbes<E> clone() {
      ListOfStatProbes<E> sa;
      try {
         sa = (ListOfStatProbes<E>)super.clone();
      }
      catch (CloneNotSupportedException cne) {
         throw new IllegalStateException ("CloneNotSupportedException for a class implementing Cloneable");
      }
      if (probes != null)
         sa.probes = new ArrayList<E> (probes);
      return sa;
   }


   public boolean add (E o) {
      return probes.add (o);
   }

   public void add (int index, E o) {
      probes.add (index, o);
   }

   public boolean addAll (Collection<? extends E> c) {
      return probes.addAll (c);
   }

   public boolean addAll (int index, Collection<? extends E> c) {
      return probes.addAll (index, c);
   }

   public void clear() {
      probes.clear();
   }

   public boolean contains (Object o) {
      return probes.contains (o);
   }

   public boolean containsAll (Collection<?> c) {
      return probes.containsAll (c);
   }

   public boolean equals (Object o) {
      return probes.equals (o);
   }

   public E get (int index) {
      return probes.get (index);
   }

   public int hashCode() {
      return probes.hashCode();
   }

   public int indexOf (Object o) {
      return probes.indexOf (o);
   }

   public boolean isEmpty() {
      return probes.isEmpty();
   }

   public Iterator<E> iterator() {
      return probes.iterator();
   }

   public int lastIndexOf (Object o) {
      return probes.lastIndexOf (o);
   }

   public ListIterator<E> listIterator() {
      return probes.listIterator();
   }

   public ListIterator<E> listIterator (int index) {
      return probes.listIterator();
   }

   public E remove (int index) {
      return probes.remove (index);
   }

   public boolean remove (Object o) {
      return probes.remove (o);
   }

   public boolean removeAll (Collection<?> c) {
      return probes.removeAll (c);
   }

   public boolean retainAll (Collection<?> c) {
      return probes.retainAll (c);
   }

   public E set (int index, E element) {
      return probes.set (index, element);
   }

   public int size() {
      return probes.size();
   }

   public List<E> subList (int fromIndex, int toIndex) {
      return probes.subList (fromIndex, toIndex);
   }

   public Object[] toArray() {
      return probes.toArray();
   }

   public <T> T[] toArray (T[] a) {
      return probes.toArray (a);
   }

}