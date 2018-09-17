/*
 * Class:        RedblackTree
 * Description:  implementation of class EventList using a red-black tree
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

import java.util.TreeMap;
import java.util.Comparator;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.LinkedList;
import java.util.ConcurrentModificationException;
import java.util.NoSuchElementException;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.simevents.Event;

/**
 * An implementation of  @ref EventList using a *red black* tree, which is
 * similar to a binary search tree except that every node is colored red or
 * black. When modifying the structure, the tree is reorganized for the
 * colors to satisfy rules that give an average @f$O(\log(n))@f$ time for
 * removing the first event or inserting a new event, where @f$n@f$ is the
 * number of elements in the structure. However, adding or removing events
 * imply reorganizing the tree and requires more overhead than a binary
 * search tree.
 *
 * The present implementation uses the Java 2  TreeMap class which implements
 * a red black tree for general usage. This event list implementation is not
 * efficient.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RedblackTree implements EventList {
   private final TreeMap<Event, Node> tree = new TreeMap<Event, Node>(new EventComparator());
   private int modCount = 0;
   
   @Override
   public void clear() {
      Iterator<Node> itr = tree.values().iterator();
      while (itr.hasNext()) {
         Node node = itr.next();
         node.events.clear();
         itr.remove();
         node.nextNode = null;
      }
      ++modCount;
   }

   @Override
   public void add (Event ev) {
      Node node = tree.get (ev);
      if (node != null)
         node.events.add (ev);
      else
         tree.put (new EventMapKey (ev), newNode (ev));
      ++modCount;
   }

   @Override
   public void addFirst (Event ev) {
      Node node = tree.get (ev);
      if (node != null)
         node.events.add (ev);
      else
         tree.put (new EventMapKey(ev), newNode (ev));
      ++modCount;
   }

   @Override
   public void addBefore (Event ev, Event other) {
      Node node = tree.get (other);
      if (node == null)
         throw new IllegalArgumentException ("Event not in list.");
      node.addBefore (ev, other);
      ++modCount;
   }

   @Override
   public void addAfter (Event ev, Event other) {
      Node node = tree.get (other);
      if (node == null)
         throw new IllegalArgumentException ("Event not in list.");
      node.addAfter (ev, other);
      ++modCount;
   }

   @Override
   public Event getFirst() {
      return isEmpty() ? null :
         tree.get (tree.firstKey()).events.get (0);
   }

   @Override
   public Event getFirstOfClass (String cl) {
      Iterator<Node> itr = tree.values().iterator();
      while (itr.hasNext()) {
         Node node = itr.next();
         Event ev = node.getFirstOfClass (cl);
         if (ev != null)
            return ev;
      }
      return null;
   }

   @Override
   public <E extends Event> E getFirstOfClass (Class<E> cl) {
      Iterator<Node> itr = tree.values().iterator();
      while (itr.hasNext()) {
         Node node = itr.next();
         E ev = node.getFirstOfClass (cl);
         if (ev != null)
            return ev;
      }
      return null;
   }

   @Override
   public boolean remove (Event ev) {
      Node node = tree.get (ev);
      if (node == null)
         return false;
      if (node.remove (ev)) {
         tree.remove (ev);
         node.nextNode = null;
      }
      ++modCount;
      return true;
   }

   @Override
   public Event removeFirst() {
      if (tree.isEmpty())
         return null;
      Event evKey = tree.firstKey();
      Node node = tree.get (evKey);
      Event first = node.events.get (0);
      node.events.remove (0);
      if (node.events.isEmpty()) {
         tree.remove (evKey);
         node.nextNode = null;
      }
      ++modCount;
      return first;
   }

   @Override
   public String toString() {
      StringBuilder sb = new StringBuilder ("Contents of the event list RedblackTree:");
      for (Node node : tree.values()) {
         for (Event ev : node.events)
            sb.append (PrintfFormat.NEWLINE +
                       PrintfFormat.g (12, 7, ev.time()) + ", " +
                       PrintfFormat.g (8, 4, ev.priority()) + " : " +
                       ev.toString());
      }
      return sb.toString();
   }

   @Override
   public Iterator<Event> iterator() {
      return listIterator();
   }

   @Override
   public ListIterator<Event> listIterator() {
      return new RBItr();
   }

   @Override
   public boolean isEmpty() {
      return tree.isEmpty();
   }

   private static class Node {
      // For the iterator
      public Node prevNode = null;
      public Node nextNode = null;

      public java.util.List<Event> events = new LinkedList<Event>();

      public Node (Event ev) {
         events.add (ev);
      }

      public void addAfter (Event ev, Event other) {
         ListIterator<Event> itr = events.listIterator();
         while (itr.hasNext()) {
            Event listev = itr.next();
            if (listev == other) {
               itr.add (ev);
               return;
            }
         }
         throw new IllegalArgumentException ("Event not in node.");
      }

      public void addBefore (Event ev, Event other) {
         ListIterator<Event> itr = events.listIterator();
         while (itr.hasNext()) {
            Event listev = itr.next();
            if (listev == other) {
               itr.previous();
               itr.add (ev);
               return;
            }
         }
         throw new IllegalArgumentException ("Event not in node.");
      }

      public Event getFirstOfClass (String cl) {
         Iterator<Event> itr = events.iterator();
         while (itr.hasNext()) {
            Event listev = itr.next();
            if (listev.getClass().getName().equals (cl))
               return listev;
         }
         return null;
      }

      @SuppressWarnings("unchecked")
      public <E extends Event> E getFirstOfClass (Class<E> cl) {
         Iterator<Event> itr = events.iterator();
         while (itr.hasNext()) {
            Event listev = itr.next();
            if (listev.getClass() == cl)
               return (E)listev;
         }
         return null;
      }

      /**
       * Remove an event from a node.
       * Returns true if the node becomes empty.
       */
      public boolean remove (Event ev) {
         Iterator<Event> itr = events.iterator();
         while (itr.hasNext()) {
            Event listev = itr.next();
            if (listev == ev) {
               itr.remove();
               return events.isEmpty();
            }
         }
         throw new IllegalArgumentException ("Event not in node.");
      }

      @Override
      public String toString() {
         StringBuilder sb = new StringBuilder();
         boolean first = true;
         Iterator<Event> itr = events.iterator();
         while (itr.hasNext()) {
            if (first)
               first = false;
            else
               sb.append (", ");
            sb.append (itr.next());
         }
         return sb.toString();
      }
   }

   private static class EventComparator implements Comparator<Event> {
      @Override
      public int compare (Event ev1, Event ev2) {
         return ev1.compareTo(ev2);
      }
   }

   private class RBItr implements ListIterator<Event> {
      private int expectedModCount;
      private Node prevNode;
      private Node nextNode;
      private int prevNodeIndex;
      private int nextNodeIndex;
      private int nextIndex;

      RBItr() {
         expectedModCount = modCount;
         prevNode = null;
         nextNode = tree.isEmpty() ? null :
            (Node)tree.get (tree.firstKey());
         prevNodeIndex = 0;
         nextNodeIndex = 0;
         nextIndex = 0;

         Iterator<Node> itr = tree.values().iterator();
         Node lastNode = null;
         while (itr.hasNext()) {
            Node node = itr.next();
            node.prevNode = lastNode;
            if (lastNode != null)
               lastNode.nextNode = node;
            node.nextNode = null;
            lastNode = node;
         }
      }

      @Override
      public void add(Event ev) {
         throw new UnsupportedOperationException();
      }

      @Override
      public boolean hasNext() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return nextNode != null &&
            nextNodeIndex < nextNode.events.size();
      }

      @Override
      public boolean hasPrevious() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return prevNode != null &&
            prevNodeIndex >= 0;
      }

      @Override
      public Event next() {
         if (!hasNext())
            throw new NoSuchElementException();

         ++nextIndex;
         Event ev = (Event)nextNode.events.get(nextNodeIndex);
         prevNode = nextNode;
         prevNodeIndex = nextNodeIndex;
         ++nextNodeIndex;
         if (nextNodeIndex >= nextNode.events.size()) {
            nextNode = nextNode.nextNode;
            nextNodeIndex = 0;
         }
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
         Event ev = (Event)prevNode.events.get(prevNodeIndex);
         nextNode = prevNode;
         nextNodeIndex = prevNodeIndex;
         --prevNodeIndex;
         if (prevNodeIndex < 0) {
            prevNode = prevNode.prevNode;
            if (prevNode != null)
               prevNodeIndex = prevNode.events.size() - 1;
         }
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
         throw new UnsupportedOperationException();
      }

      @Override
      public void set (Event ev) {
         throw new UnsupportedOperationException();
      }
   }

   private Node newNode (Event ev) {
      return new Node (ev);
   }

   private class EventMapKey extends Event {
      public EventMapKey(Event ev) {
        this.eventTime = ev.time();
        this.priority = ev.priority();
      }

      @Override
      public void actions() { }
   }
}