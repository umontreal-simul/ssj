/*
 * Class:        RedblackTree
 * Description:  implementation of class EventList using a red-black tree
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
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
import umontreal.ssj.simevents.Sim;

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
   private TreeMap<Event, Node> tree = new TreeMap<Event, Node>(new EventComparator());
   private static Node free = null;
   private int modCount = 0;
   public void clear() {
      Iterator<Node> itr = tree.values().iterator();
      while (itr.hasNext()) {
         Node node = itr.next();
         node.events.clear();
         itr.remove();
         synchronized (RedblackTree.class) {
            node.nextNode = free;
            free = node;
         }
      }
      ++modCount;
   }

   public void add (Event ev) {
      Node node = tree.get (ev);
      if (node != null)
         node.events.add (ev);
      else
         tree.put (new EventMapKey (ev), newNode (ev));
      ++modCount;
   }

   public void addFirst (Event ev) {
 //     ev.setTime (Sim.time());   // necessaire si eventime n'est pas deja a 0
      Node node = tree.get (ev);
      if (node != null)
         node.events.add (ev);
      else
         tree.put (new EventMapKey(ev), newNode (ev));
      ++modCount;
   }

   public void addBefore (Event ev, Event other) {
      Node node = tree.get (other);
      if (node == null)
         throw new IllegalArgumentException ("Event not in list.");
   //   ev.setTime (other.time());
      node.addBefore (ev, other);
      ++modCount;
   }

   public void addAfter (Event ev, Event other) {
      Node node = tree.get (other);
      if (node == null)
         throw new IllegalArgumentException ("Event not in list.");
   //   ev.setTime (other.time());
      node.addAfter (ev, other);
      ++modCount;
   }

   public Event getFirst() {
      return isEmpty() ? null :
         tree.get (tree.firstKey()).events.get (0);
   }

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

   public boolean remove (Event ev) {
      Node node = tree.get (ev);
      if (node == null)
         return false;
      if (node.remove (ev)) {
         tree.remove (ev);
         synchronized (RedblackTree.class) {
            node.nextNode = free; free = node;
         }
      }
      ++modCount;
      return true;
   }

   public Event removeFirst() {
      if (tree.isEmpty())
         return null;
      Event evKey = tree.firstKey();
      Node node = tree.get (evKey);
      Event first = node.events.get (0);
      node.events.remove (0);
      if (node.events.isEmpty()) {
         tree.remove (evKey);
         synchronized (RedblackTree.class) {
            node.nextNode = free; free = node;
         }
      }
      ++modCount;
      return first;
   }

   public String toString() {
      StringBuffer sb = new StringBuffer ("Contents of the event list RedblackTree:");
      for (Node node : tree.values()) {
         for (Event ev : node.events)
            sb.append (PrintfFormat.NEWLINE +
                       PrintfFormat.g (12, 7, ev.time()) + ", " +
                       PrintfFormat.g (8, 4, ev.priority()) + " : " +
                       ev.toString());
      }
      return sb.toString();
   }

   public Iterator<Event> iterator() {
      return listIterator();
   }

   public ListIterator<Event> listIterator() {
      return new RBItr();
   }

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

      public String toString() {
         StringBuffer sb = new StringBuffer();
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
      public int compare (Event ev1, Event ev2) {
         return ev1.compareTo(ev2);
      }

      public boolean equals (Object obj) { return true; }
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

      public void add(Event ev) {
         throw new UnsupportedOperationException();
      }

      public boolean hasNext() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return nextNode != null &&
            nextNodeIndex < nextNode.events.size();
      }

      public boolean hasPrevious() {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return prevNode != null &&
            prevNodeIndex >= 0;
      }

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

      public int nextIndex() {
         if (!hasNext())
            throw new NoSuchElementException();

         return nextIndex;
      }

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

      public int previousIndex() {
         if (!hasPrevious())
            throw new NoSuchElementException();

         return nextIndex - 1;
      }

      public void remove() {
         throw new UnsupportedOperationException();
      }

      public void set (Event ev) {
         throw new UnsupportedOperationException();
      }
   }

   private Node newNode (Event ev) {
      Node temp;
      synchronized (RedblackTree.class) {
         if (free == null)
            return new Node (ev);

         temp = free;
         free = free.nextNode;
      }
      temp.events.add (ev);
      return temp;
   }

   private class EventMapKey extends Event {
      public EventMapKey(Event ev) {
        this.eventTime = ev.time();
        this.priority = ev.priority();
      }

      public void actions() { }
   }
}