/*
 * Class:        BinaryTree
 * Description:  implementation of class EventList using a binary search tree
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
import umontreal.ssj.simevents.Event;
import umontreal.ssj.util.PrintfFormat;


/**
 * An implementation of  @ref EventList using a binary search tree. Every
 * event is stored into a tree node which has left and right children. Using
 * the event time as a comparator the left child is always smaller than its
 * parent whereas the right is greater or equal. This allows an average
 * @f$O(\log(n))@f$ time for adding an event and searching the first event,
 * where @f$n@f$ is the number of events in the structure. There is less
 * overhead for adding and removing events than splay tree or red black tree.
 * However, in the worst case, adding or removing could be done in time
 * proportional to @f$n@f$ because the binary search tree can be turned into
 * a linked list.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BinaryTree implements EventList {
   // racine de l'arbre
   private Entry root = null;

   // compteur de modifications sur l'iterateur.
   private int modCount = 0;
   
   @Override
   public boolean isEmpty() {
      return root == null;
   }

   @Override
   public void clear() {
      while (root != null)
         remove (root);
   }

   @Override
   public void add (Event ev) {
      // fonction qui ajoute un evenement dans l'arbre
      // note : si deux evenements ont le meme temps, alors il faut
      //        toujours faire en sorte que ces evenements se retrouvent
      //        comme les fils droits les uns des autres

      Entry cursor = root;
      boolean found = false;

      if (cursor == null)
         root = add (ev, null);
      else {
         while (!found) {
            if (ev.compareTo(cursor.event) < 0) {
               if (cursor.left == null) {
                  cursor.left = add (ev, cursor);
                  found = true;
               }
               cursor = cursor.left;

            } else {
               if (cursor.right == null) {
                  cursor.right = add (ev, cursor);
                  found = true;
               }
               cursor = cursor.right;
            }
         }
      }
      ++modCount;
   }

   @Override
   public void addFirst (Event ev) {
   /**
    * Ajoute "ev" comme premier evenement dans l'arbre.
    * Donc completement a gauche
    * On met l'ancien premier evenement a droite de ev.
    * (Necessaire quand on a des evenements simultanes)
    */
      Entry cursor = root;

      if (cursor != null) {
         while (cursor.left != null)
            cursor = cursor.left;

         Entry e = add (ev, cursor.father);
         e.right = cursor;
         if (cursor == root)
            root = e;
         else
            cursor.father.left = e;
         cursor.father = e;
      }
      else
         root = add (ev, null);

      ++modCount;
   }

   @Override
   public void addBefore (Event ev, Event other) {
      Entry otherEntry = findEntry (other);
      Entry evEntry    = add (ev , null);

      if (otherEntry == null)
         throw new IllegalArgumentException("other not in the tree");

      // insere evEntry a la place de otherEntry et otherEntry
      // devient le fils droit de evEntry
      if (otherEntry != root) {
         if (otherEntry == otherEntry.father.right)
            otherEntry.father.right = evEntry;
         else
            otherEntry.father.left = evEntry;
      }
      else
         root = evEntry;

      evEntry.father = otherEntry.father;
      otherEntry.father   = evEntry;
      evEntry.right   = otherEntry;

      // le ss-arbre de droite de otherEntry devient le
      // ss-arbre de droite de evEntry
      // permet que evEntry soit exactement apres
      // otherEntry qu'importe les operations effectuees
      evEntry.left  = otherEntry.left;

      if (evEntry.left != null)
         evEntry.left.father  = evEntry;

      otherEntry.left = null;

      ++modCount;
   }

   @Override
   public void addAfter (Event ev, Event other) {
      // on va chercher le "Entry" de other
      Entry otherEntry = findEntry (other);

      if (otherEntry == null)
         throw new IllegalArgumentException("other not in the tree");

      // otherEntry est le parent de evEntry
      Entry evEntry = add (ev, otherEntry);

      evEntry.right = otherEntry.right;
      otherEntry.right = evEntry;

      if (evEntry.right != null)
         evEntry.right.father = evEntry;

      ++modCount;
   }

   @Override
   public Event getFirst() {
      if (root==null)
         return null;
      Entry cursor = root;
      while (cursor.left != null)
         cursor = cursor.left;
      return cursor.event;
   }

   @Override
   public Event getFirstOfClass (String cl) {
      Entry cursor = root;
      if (root != null)
         while (cursor.left != null)
            cursor = cursor.left;

      while (cursor != null) {
         if (cursor.event.getClass().getName().equals (cl))
            return cursor.event;
         cursor = successor (cursor);
      }
      return null;
   }

   @SuppressWarnings("unchecked")
   @Override
   public <E extends Event> E getFirstOfClass (Class<E> cl) {
      Entry cursor = root;
      if (root != null)
         while (cursor.left != null)
            cursor = cursor.left;

      while (cursor != null) {
         if (cursor.event.getClass() == cl)
            return (E)cursor.event;
         cursor = successor (cursor);
      }
      return null;
   }

   @Override
   public Iterator<Event> iterator() {
      return listIterator();
   }

   @Override
   public ListIterator<Event> listIterator() {
      return new BTItr();
   }

   @Override
   public boolean remove (Event ev) {
      Entry evEntry = findEntry(ev);
      if (evEntry == null)
         return false;
      else
         return remove(evEntry);
   }

   @Override
   public Event removeFirst() {
      if (root == null)
         return null;

      Entry cursor = root;
      while (cursor.left != null)
         cursor = cursor.left;

      Event first = cursor.event;
      remove(cursor);

      return first;
   }

   @Override
   public String toString() {
      StringBuilder sb = new StringBuilder ("Contents of the event list BinaryTree:");
      Entry cursor = root;

      if (root != null)
         while (cursor.left != null)
            cursor = cursor.left;

      while (cursor != null) {
         sb.append (PrintfFormat.NEWLINE +
                    PrintfFormat.g (12, 7, cursor.event.time()) + ", " +
                    PrintfFormat.g (8, 4, cursor.event.priority()) +
                    " : " + cursor.event.toString());
         cursor = successor (cursor);
      }

      return sb.toString();
   }

   private Entry add (Event ev, Entry father) {
      return new Entry(ev, null, null, father);
   }

   private boolean remove (Entry e) {
      boolean filsGauche = false;
      boolean isRoot = false;
      Entry cursor;

      if (e == root)
         isRoot = true;
      else {
         if (e == e.father.left)
            filsGauche = true;
         else
            filsGauche = false;
      }

      // Si condition vrai, a un fils droit ou rien
      if (e.left == null) {
         if (isRoot)
            root = e.right;
         else if (filsGauche)
            e.father.left = e.right;
         else
            e.father.right = e.right;

         if (e.right != null)
            e.right.father = e.father;
      }
      else if (e.right == null) {
         // Si condition vrai,  a uniquement un fils gauche
         if (isRoot)
            root = e.left;
         else if (filsGauche)
            e.father.left = e.left;
         else
            e.father.right = e.left;
         e.left.father = e.father;
      }
      else {
         // a 2 fils
         // recherche son descendant le plus petit dans le ss-arbre de droite
         // et remplace "e" par ce descendant

         cursor = e.right;
         if (cursor.left == null) {
            // c'est son fils de droite
            if (isRoot)
               root = cursor;
            else {
               if (filsGauche)
                  e.father.left = cursor;
               else
                  e.father.right = cursor;
            }
            cursor.left = e.left;
         }
         else {
            // recherche de la plus petite valeur dans le ss-arbre droit
            while (cursor.left != null)
               cursor = cursor.left;

            if (isRoot)
               root = cursor;
            else if (filsGauche)
               e.father.left = cursor;
            else
               e.father.right = cursor;

            // echange entre e et cursor et elimination de e
            cursor.father.left = cursor.right;
            if (cursor.right != null)
               cursor.right.father = cursor.father;

            cursor.right = e.right;
            cursor.left  = e.left;
            e.right.father = cursor;
         }

         cursor.father = e.father;
         e.left.father = cursor;
      }

      // recupere l'espace du noeud
      e.right = null;
      e.left =  null;
      e.event = null;

      ++modCount;
      return true;
   }

   private Entry successor (Entry cursor) {
      if (cursor == null)
         return null;

      if (cursor.right != null) {
         cursor = cursor.right;
         while (cursor.left != null)
            cursor = cursor.left;
      }
      else {
         while (cursor.father != null && cursor.father.right == cursor)
            cursor = cursor.father;
         cursor = cursor.father;
      }
      return cursor;
   }

   /**
    * fonction qui trouve le noeud (Entry) d'un evenement
    * dans l'arbre
    */
   private Entry findEntry (Event ev) {
      Entry cursor = root;
      while (cursor != null) {
         if (cursor.event == ev)
            return cursor;
         else if (ev.compareTo(cursor.event) < 0)
            cursor = cursor.left;
         else
            cursor = cursor.right;
      }
      return null;
   }

   private Entry predecessor (Entry cursor) {
      if (cursor == null)
         return null;

      if (cursor.left != null) {
         cursor = cursor.left;
         while (cursor.right != null)
            cursor = cursor.right;
      }
      else {
         while (cursor.father != null && cursor.father.left == cursor)
            cursor = cursor.father;
         cursor = cursor.father;
      }
      return cursor;
   }

   /**
    * Classe interne representant les noeuds de l'arbre
    */
   private static class Entry  {
      Event event;
      Entry right;
      Entry left;
      Entry father;

      Entry (Event event, Entry left, Entry right, Entry father) {
         this.event = event;
         this.left = left;
         this.right = right;
         this.father = father;
      }
   }

   private class BTItr implements ListIterator<Event> {
      private Entry prev;
      private Entry next;
      private Entry lastRet;
      private int expectedModCount;
      private int nextIndex;

      BTItr() {
         prev = null;
         next = root;
         if (next != null) {
            while (next.left != null)
               next = next.left;
         }
         expectedModCount = modCount;
         lastRet = null;
         nextIndex = 0;
      }

      @Override
      public void add(Event ev) {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();

         // Check the event time and priority
         if (next != null && ev.compareTo(next.event) > 0) {
            ev.setTime (next.event.time());
            ev.setPriority (next.event.priority());
         }
         if (prev != null && ev.compareTo(prev.event) < 0) {
            ev.setTime (prev.event.time());
            ev.setPriority (prev.event.priority());
         }

         Entry e = BinaryTree.this.add (ev, next);
         if (prev != null) {
            // Ajouter ev apr`es prev.
            // insere e comme fils droit de prev
            e.father = prev;
            e.right = prev.right;
            prev.right = e;
            if (e.right != null)
               e.right.father = e;
         }
         else {
            // ajoute ev avant next.
            // insere e a la place de eo et eo devient le fils droit de e
            if (next != root) {
               if (next == next.father.left)
                  next.father.left = e;
               else
                  next.father.right = e;
            }
            else
               root = e;
            e.father = prev.father;
            prev.father = e;
            e.left = prev;
            // le ss-arbre de droite de eo devient le ss-arbre de droite de e
            // permet que e soit exactement apres eo qu'importe les
            // operations effectuees
            e.right = prev.right;
            if (e.right != null)
               e.right.father = e;
            prev.right = null;
         }

         prev = e;
         ++nextIndex;
         lastRet = null;
         ++modCount;
         ++expectedModCount;
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
         Event ev = next.event;
         lastRet = next;
         prev = next;
         next = successor (next);
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
         Event ev = prev.event;
         lastRet = prev;
         next = prev;
         prev = predecessor (prev);
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
            next = successor (next);
         else { // Last call to next or no call
            prev = predecessor (prev);
            --nextIndex;
         }
         BinaryTree.this.remove (lastRet);
         lastRet = null;
         ++expectedModCount;
      }

      @Override
      public void set (Event ev) {
         if (modCount != expectedModCount)
            throw new ConcurrentModificationException();
         if (lastRet == null)
            throw new IllegalStateException();

         Entry pred = predecessor (lastRet);
         Entry succ = successor (lastRet);
         if (pred != null && ev.compareTo(pred.event) < 0) {
            ev.setTime (pred.event.time());
            ev.setPriority (pred.event.priority());
         }
         if (succ != null && ev.compareTo(succ.event) > 0) {
            ev.setTime (succ.event.time());
            ev.setPriority (succ.event.priority());
         }
         lastRet.event = ev;
      }
   }

/*   public static void main (String[] args) {
      BinaryTree sp = new BinaryTree();

      Event1 e1 = new Event1(); e1.setTime(10.0);
      Event1 e2 = new Event1(); e2.setTime(20.0);
      Event1 e3 = new Event1(); e3.setTime(30.0);
      Event1 e4 = new Event1(); e4.setTime(40.0);
      Event1 e5 = new Event1(); e5.setTime(50.0);
      Event1 e6 = new Event1(); e6.setTime(60.0);
      Event1 e7 = new Event1(); e7.setTime(70.0);

      sp.add(e1);
      sp.add(e2);
      sp.print2(sp.root);
      sp.add(e3);
      sp.print2(sp.root);
      sp.add(e4);
      sp.print2(sp.root);
      sp.add(e5);
      sp.print2(sp.root);
      sp.add(e6);
      sp.print2(sp.root);
      sp.add(e7);
      sp.print2(sp.root);
      // sp.add(e5);
      // sp.print2(sp.root);
      sp.add(e7);
      sp.print2(sp.root);
      // sp.add(e5);
      // sp.print2(sp.root);

      sp.getFirst();
      System.out.println(".....after GetFirst" +
                         PrintfFormat.NEWLINE +
                         PrintfFormat.NEWLINE +
                         PrintfFormat.NEWLINE);
      sp.print2(sp.root);
      sp.remove(e3);
      System.out.println("Apres remove" + PrintfFormat.NEWLINE);
      sp.print2(sp.root);
   }

   private void print(Entry t) {
      if (t != null){
         print (t.left);
         System.out.println ("===========> Event time "+t.event.time());
         print (t.right);
      }
   }

   private void print2(Entry t) {
      System.out.println("===============================  "+
                         "print2 : pour ..... "+t.event.time());
      if (t != null) {
         System.out.println ("===========> ev time   "+t.event.time());
         gauche (t.left);
         droite (t.right);
         System.out.println();
      }
      else
         System.out.println ("===========> gauche  null ");
   }

   private void gauche (Entry t) {
      if (t != null) {
         System.out.println ("===========> gauche   "+t.event.time());
         gauche (t.left);
         droite (t.right);
         System.out.println();
      }
      else
         System.out.println ("===========> gauche  null ");
   }

   private void droite (Entry t) {
      if (t != null){
         System.out.println ("===========> droite  "+t.event.time());
         gauche (t.left);
         droite (t.right);
         // System.out.println();
      }
      else
         System.out.println ("===========> droite  null ");
   }

   private static class Event1 extends Event {
      public void actions() {}

      public String toString() {
         return "Event(" + eventTime + ")";
      }
   };





*/
   /*
        public  BinaryTree() {

      Event1 e1 = new Event1(); e1.setTime(7.0);
      Event1 e2 = new Event1(); e2.setTime(5.0);
      Event1 e3 = new Event1(); e3.setPriority(2); e3.setTime(5.0);
      Event1 e4 = new Event1(); e4.setTime(6.0);
      Event1 e5 = new Event1(); e5.setPriority(2); e5.setTime(10.0);
      Event1 e6 = new Event1(); e6.setTime(9.0);
      Event1 e7 = new Event1(); e7.setTime(11.0);

        add(e1);
        add(e2);
        add(e3);
        add(e4);
        add(e5);
        add(e6);
        add(e7);
        print22(root);
        remove (e5);
        print22(root);

        }



        public void print22(Entry t) {
        System.out.println("racine............ ..... "+t.event.time());
        gauche2(t.left);
        droite2(t.right);
        System.out.println();

        }



        public void gauche2 (Entry t) {
        if (t!=null){
        System.out.println ("===========> gauche   "+t.event.time());
        gauche2(t.left);
        droite2(t.right);
        System.out.println();

        }
        else System.out.println ("===========> gauche  null ");
        }

        public void droite2 (Entry t) {
        if (t!=null){
        System.out.println ("===========> droite  "+t.event.time());
        gauche2(t.left);
        droite2(t.right);
        // System.out.println();

        }
        else System.out.println ("===========> droite  null ");
        }
   */
}