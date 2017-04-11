/*
 * Class:        DoublyLinked
 * Description:  implementation of class EventList using a doubly-linked list
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
package umontreal.ssj.simevents.eventlist;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.ConcurrentModificationException;
import java.util.NoSuchElementException;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.simevents.Event;

/**
 * An implementation of  @ref EventList using a doubly linked linear list.
 * Each event is stored into a list node that contains a pointer to its
 * following and preceding events. Adding an event requires a linear search
 * to keep the event list sorted by event time and priority. Removing the
 * first event is done in constant time because it simply removes the first
 * list node. List nodes are recycled for increased memory management
 * efficiency.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class DoublyLinked implements EventList {
   private int modCount = 0;

   // First and last elements in the list.
   private Node first = null, last = null;

   @Override
   public boolean isEmpty() {
      return first == null;
   }

   @Override
   public void clear() {
      if (first == null)
         return;

      last.succ = null;
      last = first = null;
      ++modCount;
   }

   @Override
   public void add (Event ev) {
      Node newNode;
      newNode = new Node();

      newNode.ev = ev;
      ++modCount;
      if (last == null) {     // Easy: the event list was empty.
         first = last = newNode;
         first.prec = first.succ = null;
         return;
      }

      Node node = findPosition (ev);  // Finds where to insert.
      if (node == null) {             // Must be inserted first.
         newNode.succ = first;
         newNode.succ.prec = newNode;
         first = newNode;
         newNode.prec = null;
      }
      else {                          // Insert after node.
         newNode.prec = node;
         newNode.succ = node.succ;
         node.succ = newNode;
         if (newNode.succ != null)
            newNode.succ.prec = newNode;
         else
            last = newNode;
      }
   }

   @Override
   public void addFirst (Event ev) {
      Node newNode;
      newNode = new Node();
      
      newNode.ev = ev;
      newNode.prec = null;
      if (first == null) {
         first = last = newNode;
         first.succ = null;
      }
      else {
         newNode.succ = first;
         first.prec = newNode;
         first = newNode;
      }
      ++modCount;
   }

   @Override
   public void addBefore (Event ev, Event other) {
      Node node = last;
      while (node != null && node.ev.compareTo(other) >= 0 && node.ev != other)
         node = node.prec;
      if (node.ev != other)
         throw new IllegalArgumentException ("Event not in list.");

      Node newNode;
      newNode = new Node();
      
      newNode.ev = ev;

      newNode.prec = node.prec;
      newNode.succ = node;
      node.prec = newNode;
      if (newNode.prec != null)
         newNode.prec.succ = newNode;
      else
         first = newNode;
      ++modCount;
   }

   @Override
   public void addAfter (Event ev, Event other) {
      Node node = last;
      while (node != null && node.ev.compareTo(other) >= 0 && node.ev != other)
         node = node.prec;
      if (node.ev != other)
         throw new IllegalArgumentException ("Event not in list.");

      Node newNode;
      newNode = new Node();
      
      newNode.ev = ev;

      newNode.prec = node;
      newNode.succ = node.succ;
      node.succ = newNode;
      if (newNode.succ != null)
         newNode.succ.prec = newNode;
      else
         last = newNode;
      ++modCount;
   }

   @Override
   public Event getFirst() {
      return first == null ? null : first.ev;
   }

   @Override
   public Event getFirstOfClass (String cl) {
      Node node = first;
      while (node != null) {
         if (node.ev.getClass().getName().equals (cl))
            return node.ev;
         node = node.succ;
      }
      return null;
   }

   @SuppressWarnings("unchecked")
   @Override
   public <E extends Event> E getFirstOfClass (Class<E> cl) {
      Node node = first;
      while (node != null) {
         if (node.ev.getClass() == cl)
            return (E)node.ev;
         node = node.succ;
      }
      return null;
   }

   @Override
   public boolean remove (Event ev) {
      // Find the node corresponding to this event ev.
      Node node = last;
      while (node != null && node.ev.compareTo(ev) >= 0 && node.ev != ev)
         node = node.prec;
      if (node == null || node.ev != ev)
         return false;

      if (node == last && node == first)
         last = first = null; // The list is now empty.
      else {
         if (node == last) {
            last = node.prec;
            last.succ = null;
         }
         else
            node.succ.prec = node.prec;
         if (node == first) {
            first = node.succ;
            first.prec = null;
         }
         else {
            node.prec.succ = node.succ;
            node.prec = null;
         }
      }
      node.ev = null;
      node.succ = null;
              
      ++modCount;
      return true;
   }

   @Override
   public Event removeFirst() {
      if (first == null)
         return null;

      Event ev = first.ev;
      Node temp = first;
      first = first.succ;
      if (first == null)
         last = null;
      else
         first.prec = null;

      temp.ev = null;
      temp.succ = null;

      ++modCount;
      return ev;
   }

   @Override
   public Iterator<Event> iterator() {
      return listIterator();
   }

   @Override
   public ListIterator<Event> listIterator() {
      return new DLItr();
   }

   @Override
   public String toString() {
      StringBuilder sb = new StringBuilder ("Contents of the event list DoublyLinked:");
      Node node = first;
      while (node != null) {
         sb.append (PrintfFormat.NEWLINE +
                    PrintfFormat.g (12, 7, node.ev.time()) + ", " +
                    PrintfFormat.g (8, 4, node.ev.priority()) + " : " +
                    node.ev.toString());
         node = node.succ;
      }
      return sb.toString();
   }

   // A element of the event list. This node contains the event ev.
   // His predecessor and successor are prec and succ.
   private static class Node {
      Event ev;
      Node prec, succ;
   }

   private class DLItr implements ListIterator<Event> {
      private Node prev;
      private Node next;
      private Node lastRet;
      private int expectedModCount;
      private int nextIndex;

      DLItr() {
         prev = null;
         next = first;
         expectedModCount = modCount;
         lastRet = null;
         nextIndex = 0;
      }

      @Override
      public void add(Event ev) {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();

         // Check the event time and priority
         if (next != null && ev.compareTo(next.ev) > 0) {
            ev.setTime (next.ev.time());
            ev.setPriority (next.ev.priority());
         }
         if (prev != null && ev.compareTo(prev.ev) < 0) {
            ev.setTime (prev.ev.time());
            ev.setPriority (prev.ev.priority());
         }

         Node newNode;
         newNode = new Node();
         
         newNode.ev = ev;
         ++nextIndex;
         ++modCount;
         ++expectedModCount;
         lastRet = null;
         if (last == null) {     // Easy: the event list was empty.
            first = last = newNode;
            first.prec = first.succ = null;
            prev = newNode;
            next = null;
            nextIndex = 1;
         }
         else if (prev == null) {             // Must be inserted first.
            // next is non-null or the list would be empty.
            newNode.succ = first;
            newNode.succ.prec = newNode;
            first = newNode;
            newNode.prec = null;
            prev = newNode;
         }
         else {                          // Insert after node.
            // prev is non-null but next can be null.
            newNode.prec = prev;
            newNode.succ = next;
            prev.succ = newNode;
            if (newNode.succ != null)
               newNode.succ.prec = newNode;
            else
               last = newNode;
            prev = newNode;
         }
      }

      @Override
      public boolean hasNext() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return next != null;
      }

      @Override
      public boolean hasPrevious() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return prev != null;
      }

      @Override
      public Event next() {
         if (!hasNext())
            throw new NoSuchElementException();

         ++nextIndex;
         Event ev = next.ev;
         lastRet = next;
         prev = next;
         next = next.succ;
         return ev;
      }

      @Override
      public int nextIndex() {
         if (!hasNext())
            throw new NoSuchElementException();

         return nextIndex;
      }

      @Override
      public Event previous() {
         if (!hasPrevious())
            throw new NoSuchElementException();

         --nextIndex;
         Event ev = prev.ev;
         lastRet = prev;
         next = prev;
         prev = prev.prec;
         return ev;
      }

      @Override
      public int previousIndex() {
         if (!hasPrevious())
            throw new NoSuchElementException();

         return nextIndex - 1;
      }

      @Override
      public void remove() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         if (lastRet == null)
            throw new IllegalStateException();

         if (lastRet == next) // Last call to previous
            next = next.succ;
         else { // Last call to next
            prev = prev.prec;
            --nextIndex;
         }
         if (lastRet == last && lastRet == first) {
            last = first = null; // The list is now empty.
            next = prev = null;
         }
         else {
            if (lastRet == last) {
               last = lastRet.prec;
               last.succ = null;
            }
            else
               lastRet.succ.prec = lastRet.prec;
            if (lastRet == first) {
               first = lastRet.succ;
               first.prec = null;
            }
            else {
               lastRet.prec.succ = lastRet.succ;
               lastRet.prec = null;
            }
         }
         lastRet.ev = null;
         lastRet.succ = null;

         lastRet = null;
         ++modCount;
         ++expectedModCount;
      }

      @Override
      public void set (Event ev) {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         if (lastRet == null)
            throw new IllegalStateException();

         // Check the event time and priority
         if (lastRet.prec != null && ev.compareTo(lastRet.prec.ev) < 0) {
            ev.setTime (lastRet.prec.ev.time());
            ev.setPriority (lastRet.prec.ev.priority());
         }
         if (lastRet.succ != null && ev.compareTo(lastRet.succ.ev) > 0) {
            ev.setTime (lastRet.succ.ev.time());
            ev.setPriority (lastRet.succ.ev.priority());
         }

         lastRet.ev = ev;
      }
   }

   private Node findPosition (Event ev)  {
      // This implementation appears quite inefficient  !!!!
      // Must try to improve.
      Node node = last;

      // Finds the occurrence time of the new event (evTime).
      while (node != null && ev.compareTo(node.ev) < 0) {
         node = node.prec;
      }
      return node;
   }
}