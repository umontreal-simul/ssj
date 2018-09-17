/*
 * Class:        TransformingList
 * Description:  List that dynamically transforms the elements of another list.
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
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
package umontreal.ssj.util;

import java.util.AbstractList;
import java.util.List;
import java.util.Iterator;
import java.util.ListIterator;

/**
 * Represents a list that dynamically transforms the elements of another
 * list. This abstract class defines a list containing an inner list of
 * elements of a certain type, and provides facilities to convert these inner
 * elements to outer elements of another type. A concrete subclass simply
 * needs to provide methods for converting between the inner and the outer
 * types.
 *  @param <IE>         the inner type of the elements.
 *  @param <OE>         the type of the outer elements.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class TransformingList<OE,IE> extends AbstractList<OE> {
   private List<IE> fromList;

/**
 * Creates a new transforming list wrapping the inner list `fromList`.
 *  @param fromList     the inner list.
 */
public TransformingList (List<IE> fromList) {
      this.fromList = fromList;
   }


   public List<IE> getInnerList() {
      return fromList;
   }

/**
 * Converts an element in the inner list to an element of the outer type.
 *  @param e            the inner element.
 *  @return the outer element.
 */
public abstract OE convertFromInnerType (IE e);

   /**
    * Converts an element of the outer type to an element for the inner
    * list.
    *  @param e            the outer element.
    *  @return the inner element.
    */
   public abstract IE convertToInnerType (OE e);


   public void add (int index, OE element) {
      IE fe = convertToInnerType (element);
      fromList.add (fe);
   }

   @Override
   public void clear () {
      fromList.clear();
   }

   @Override
   public OE get (int index) {
      return convertFromInnerType (fromList.get (index));
   }
   
   public Iterator<OE> iterator() {
      return new MyIterator (fromList.iterator());
   }

   public ListIterator<OE> listIterator() {
      return new MyListIterator (fromList.listIterator());
   }

   public ListIterator<OE> listIterator (int index) {
      return new MyListIterator (fromList.listIterator (index));
   }

   @Override
   public OE remove (int index) {
      return convertFromInnerType (fromList.remove (index));
   }

   @Override
   public OE set (int index, OE element) {
      IE from = convertToInnerType (element);
      from = fromList.set (index, from);
      return convertFromInnerType (from);
   }

   @Override
   public int size () {
      return fromList.size();
   }

   private class MyIterator implements Iterator<OE> {
      private Iterator<IE> itr;

      public MyIterator (Iterator<IE> itr) {
         this.itr = itr;
      }

      public boolean hasNext() {
         return itr.hasNext();
      }

      public OE next() {
         return convertFromInnerType (itr.next());
      }
      
      public void remove() {
         itr.remove();
      }
   }

   private class MyListIterator implements ListIterator<OE> {
      private ListIterator<IE> itr;
   
      public MyListIterator (ListIterator<IE> itr) {
         this.itr = itr;
      }

      public void add (OE o) {
         IE fe = convertToInnerType (o);
         itr.add (fe);
      }

      public boolean hasNext() {
         return itr.hasNext();
      }

      public boolean hasPrevious() {
         return itr.hasPrevious();
      }

      public OE next() {
         return convertFromInnerType (itr.next());
      }

      public int nextIndex() {
         return itr.nextIndex();
      }

      public OE previous() {
         return convertFromInnerType (itr.previous());
      }

      public int previousIndex() {
         return itr.previousIndex();
      }

      public void remove() {
         itr.remove();
      }

      public void set (OE o) {
         IE fe = convertToInnerType (o);
         itr.set (fe);
      }
   }

}