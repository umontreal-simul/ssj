/*
 * Class:        SplayTree
 * Description:  implementation of class EventList using a splay tree 
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
 * An implementation of  @ref EventList using a splay tree
 * @cite iSLE85a&thinsp;. This tree is like a binary search tree except that
 * when it is modified, the affected node is moved to the top. The
 * rebalancing scheme is simpler than for a *red black* tree and can avoid
 * the worst case of the linked list. This gives a @f$O(\log(n))@f$ average
 * time for adding or removing an event, where @f$n@f$ is the size of the
 * event list.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class SplayTree implements EventList {
   private Entry root = null;
   private int modCount = 0;

   private int myCompareTo (Event ev, Event other) {
      // A new event must always occur after those with the same time and
      // same priority in the Event list. getRa is used to ensure that.
      int j = ev.compareTo(other);
      if (0 != j)
         return j;
      if (ev.getRa() < other.getRa())
         return -1;
      if (ev.getRa() > other.getRa())
         return 1;
      return 0;
   }
   @Override
   public boolean isEmpty() {
      return root == null;
   }

   @Override
   public void clear() {
      // Simply root = null would be more efficient but the
      // entries would not be recuperated.
      while (root != null)
         remove (root);
   }

   @Override
   public void add (Event ev) {
      //ajoute un element dans l'arbre en faisant un "splay"

      //reference pour le splay :
      //WEISS, Mark Allen, Data Structures & Problem Solving using JAVA
      //section 22.5 Top-Down Splay Tree

      //"zig" : rotation entre un element et son parent
      //"zig-zag" simplifie : rotation entre un element et son parent
      //"zig-zig" : rotation entre le parent et le grand-parent,
      //            suivis par une rotation entre le parent et l'element

      //NOTE : on considere toujours qu'un nouvel element doit etre place
      //       apres les evenements de meme temps et meme priorite

      ev.setRa(1);    // Temporary; make sure ev is after other events with
                      // the same time and priority. Others have ra = 0.

      Entry next = root;
      // le nouvel evenement devient la racine
      Entry e = add (ev, null);
      root = e;
      if (next != null) {
         Entry left = root;
         Entry right = root;
         Entry temp = null;
         boolean end_splay = false;
         while (!end_splay) {
            if (myCompareTo(ev, next.event) > 0) {
               temp = next.right;
               if (temp == null) {
                  //cas "zig"
                  left.right = next;
                  next.father = left;
                  right.left = null;
                  end_splay = true;
               }
               else if (myCompareTo(ev, temp.event) < 0) {
                  //cas "zig-zag" simplifie
                  left.right = next;
                  next.father = left;
                  left = next;
                  next = temp;
               }
               else {
                  //cas "zig-zig"
                  next.right = temp.left;
                  if (temp.left != null)
                     temp.left.father = next;
                  left.right = temp;
                  temp.father = left;
                  temp.left = next;
                  next.father = temp;
                  left = temp;
                  next = temp.right;
                  if (next == null) {
                     right.left = null;
                     end_splay = true;
                  }
               }

            } else {
               temp = next.left;
               if (temp == null) {
                  //cas "zig"
                  right.left = next;
                  next.father = right;
                  left.right = null;
                  end_splay = true;
               }
               else if (myCompareTo(ev, temp.event) > 0) {
                  //cas "zig-zag" simplifie
                  right.left = next;
                  next.father = right;
                  right = next;
                  next = temp;
               }
               else {
                  //cas "zig-zig"
                  next.left = temp.right;
                  if (temp.right != null)
                     temp.right.father = next;
                  right.left = temp;
                  temp.father = right;
                  temp.right = next;
                  next.father = temp;
                  right = temp;
                  next = temp.left;
                  if (next == null) {
                     left.right = null;
                     end_splay = true;
                  }
               }
            }
         }
         temp = e.left;
         e.left = e.right;
         e.right = temp;
      }
      ev.setRa(0);
      ++modCount;
   }

   @Override
   public void addFirst (Event ev) {
      if (root == null)
         root = add (ev, null);
      else {
         Entry cursor = root;
         while (cursor.left != null)
            cursor = cursor.left;
         cursor.left = add (ev, cursor);
      }
      ++modCount;
   }

   @Override
   public void addBefore (Event ev, Event other) {
      Entry otherEntry = findEntry (other);
      if (otherEntry == null)
         throw new IllegalArgumentException ("Event not in list.");

      Entry e = add (ev, null);
      // insere e a la place de otherEntry et otherEntry
      // devient le fils droit de e
      if (otherEntry != root) {
         if (otherEntry == otherEntry.father.left)
            otherEntry.father.left = e;
         else
            otherEntry.father.right = e;
      }
      else
         root = e;
      e.father = otherEntry.father;
      otherEntry.father = e;
      e.right = otherEntry;
      // le sous-arbre de droite de otherEntry devient le sous-arbre de
      // droite de e.
      // permet que e soit exactement apres otherEntry qu'importe
      // les operations effectuees
      e.left = otherEntry.left;
      if (e.left != null)
         e.left.father = e;
      otherEntry.left = null;
      ++modCount;
   }

   @Override
   public void addAfter (Event ev, Event other) {
      Entry otherEntry = findEntry (other);
      if (otherEntry == null)
         throw new IllegalArgumentException ("Event not in list.");
      Entry e = add (ev, otherEntry);
      // insere e comme fils droit de otherEntry
      e.right = otherEntry.right;
      otherEntry.right = e;
      if (e.right != null)
         e.right.father = e;
      ++modCount;
   }

   @Override
   public Event getFirst() {
      if (root == null)
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
      return new SPItr();
   }

   @Override
   public boolean remove (Event ev) {
      //on trouve le noeud correspondant a l'evenement pour l'enlever
      if (root == null)
         return false;
      Entry e = findEntry (ev);
      if (e == null)
         return false;
      return remove (e);
   }

   @Override
   public Event removeFirst() {
      if (root == null)
         return null;
      Entry cursor = root;
      while (cursor.left != null)
         cursor = cursor.left;
      Event first = cursor.event;
      remove (cursor);
      return first;
   }

   @Override
   public String toString() {
      StringBuilder sb = new StringBuilder
          ("Contents of the event list SplayTree:");
      Entry cursor = root;
      if (root != null)
         while (cursor.left != null)
            cursor = cursor.left;
      while (cursor != null) {
         sb.append (PrintfFormat.NEWLINE +
                    PrintfFormat.g (12, 7, cursor.event.time()) + ", " +
                    PrintfFormat.g (8, 4, cursor.event.priority()) + " : " +
                    cursor.event.toString());
         cursor = successor (cursor);
      }
      return sb.toString();
   }

   private static class Entry {
      Event event;
      Entry father;
      Entry left;
      Entry right;

      Entry (Event event, Entry left, Entry right, Entry father) {
         this.event = event;
         this.left = left;
         this.right = right;
         this.father = father;
      }
   }

   private class SPItr implements ListIterator<Event> {
      private Entry prev;
      private Entry next;
      private Entry lastRet;
      private int expectedModCount;
      private int nextIndex;

      SPItr() {
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

         Entry e = SplayTree.this.add (ev, null);
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
            // insere e a la place de otherEntry et
            // otherEntry devient le fils droit de e
            if (next != root) {
               if (next == next.father.left)
                  next.father.left = e;
               else
                  next.father.right = e;
            }
            else
               root = e;
            e.father = next.father;
            next.father = e;
            e.right = next;
            // le sous-arbre de droite de otherEntry
            // devient le sous-arbre de droite de e.
            // permet que e soit exactement apres otherEntry qu'importe les
            // operations effectuees
            e.left = next.left;
            if (e.left != null)
               e.left.father = e;
            next.left = null;
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
         SplayTree.this.remove (lastRet);
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

   /**
    * Creates a new entry object.
    */
   private Entry add (Event ev, Entry father) {
      return new Entry (ev, null, null, father);
   }


   /*
     Remonte a la racine une entree avec un splay de type "bottom-up",
     qui prend pour parametre le noeud a remonter.

     Fonctionnement du splay :
     Jusqu'a ce que l'element devienne la racine, on fait :
     - Si le pere de l'element est la racine, on fait la rotation
       entre le pere et l'element (cas "zig").
     - Si le pere de l'element et l'element ne sont pas des fils du
       meme cote (gauche-droite ou droite-gauche), alors on fait
       la rotation entre le pere et l'element (cas "zig-zag" simplifie).
     - Si le pere de l'element et l'element sont des fils du meme cote
       (gauche-gauche ou droite-droite), alors on fait la rotation
       entre le grand-pere et le pere avant de faire la rotation entre
       le pere et l'element (cas "zig-zig").
   */
   private void splay(Entry e)
   {
      boolean left;

      Entry f;     // = father             le pere de l'element e
      Entry gf;    // = grandFather        le grand-pere de l'element e
      Entry ggf;   // = grandGrandFather   l'arriere-grand-pere de l'element e

      while(e.father != null) {
         f = e.father;
         gf = f.father;
         left = e == f.left;

         if (left)
            // cas du fils gauche
            if (gf == null) {
               // cas "zig", on fait la rotation de f (la racine) et e
               f.father = e;
               f.left = e.right;
               if(f.left != null)
                  f.left.father = f;
               e.right = f;
               e.father = null;
            }
            else if (gf.right == f) {
               // cas "zig-zag", simplifie, pareil que le cas "zig"
               gf.right = e;

               f.father = e;
               f.left = e.right;
               if(f.left != null)
                  f.left.father = f;
               e.right = f;
               e.father = gf;
            }
            else {
               // cas "zig-zig", on fait la rotation de gf avec
               // f, suivis de la rotation de e avec f
               ggf = gf.father;

               gf.left = f.right;
               if(gf.left != null)
                  gf.left.father = gf;
               f.right = gf;
               gf.father = f;

               f.left = e.right;
               if(f.left != null)
                  f.left.father = f;
               f.father = e;
               e.right = f;

               // on rattache e a son nouveau pere
               e.father = ggf;
               if(ggf != null)
                  if(ggf.left == gf)
                     ggf.left = e;
                  else
                     ggf.right = e;
            }
         else
            //cas du fils droit
            if(gf == null) {
               // cas "zig", on fait la rotation de f (la racine) et e

               f.father = e;
               f.right = e.left;
               if(f.right != null)
                  f.right.father = f;
               e.left = f;
               e.father = null;
            }
            else if(gf.left == f) {
               // cas "zig-zag", simplifie, pareil que le cas "zig"
               gf.left = e;

               f.father = e;
               f.right = e.left;
               if(f.right != null)
                  f.right.father = f;
               e.left = f;
               e.father = gf;
            }
            else {
               // cas "zig-zig", on fait la rotation de gf avec
               // f, suivis de la rotation de e avec f
               ggf = gf.father;

               gf.right = f.left;
               if(gf.right != null)
                  gf.right.father = gf;
               f.left = gf;
               gf.father = f;

               f.right = e.left;
               if(f.right != null)
                  f.right.father = f;
               f.father = e;
               e.left = f;

               // on rattache e a son nouveau pere
               e.father = ggf;
               if(ggf != null)
                  if(ggf.left == gf)
                     ggf.left = e;
                  else
                     ggf.right = e;
            }
      }
   }


   /**
      Enleve l'entree e de l'arbre.

      Procedure :
      On fait d'abord un splay sur l'entree a enlever pour la faire monter
      a la racine. Une fois a la racine, on l'enleve. On cherche ensuite
      l'element minimal du sous-arbre droit (qui contient les elements
      qui sont plus grand que e) et, avec un splay, on le remonte a la
      racine de son sous-arbre. Cet element, puisqu'il est le minimum de
      son sous-arbre, n'a pas de fils gauche. On lui rattache donc comme
      fils gauche le sous-arbre gauche.

      La raison de cette procedure plus compliquee que l'equivalent dans
      un arbre binaire ordinaire est que, pour que la garantie de
      temps d'acces amorti en O(log n) du splay tree s'applique, il
      faut que chaque acces (qui necessite une recherche) d'un element
      s'accompagne du splay de cet element a la racine.
    */
   private boolean remove(Entry e)
   {
      if (root == null || e == null)
         return false;

      splay(e);
      // e est implicitement la racine, meme si root != e

      Entry leftTree = e.left;
      Entry rightTree = e.right;

      if(leftTree != null)
         leftTree.father = null;
      if(rightTree != null)
         rightTree.father = null;

      e.left = null;
      e.event = null;
      e.right = null;

      if(rightTree == null)
         root = leftTree;
      else if(leftTree == null)
         root = rightTree;
      else {
         // on cherche le plus petit element du sous-arbre de droite
         Entry newRoot = rightTree;
         while(newRoot.left != null)
            newRoot = newRoot.left;

         // on monte newRoot en haut du sous-arbre de droite
         splay(newRoot);

         // on rattache les deux sous-arbres
         newRoot.left = leftTree;
         leftTree.father = newRoot;
         root = newRoot;
      }

      ++modCount;

      return true;
   }


   /**
    * Finds the corresponding entry for an event.
    */
   private Entry findEntry (Event ev) {
      Entry cursor = root;
      while (cursor != null) {
         if (cursor.event == ev)
            return cursor;
         else if (ev.compareTo(cursor.event) > 0)
            cursor = cursor.right;
         else if (ev.compareTo(cursor.event) < 0)
            cursor = cursor.left;
         else {
            // on a trouve un element "egal" (compareTo() revoie 0 ), on ne peut donc plus
            // utiliser la recherche binaire et il faut donc regarder
            // tous les elements "egaux" qui se retrouvent soit apres,
            // soit avant l'element trouve
            Entry center = cursor;

            // recherche avant center
            while (cursor != null && ev.compareTo(cursor.event) == 0)
               if(cursor.event == ev)
                  return cursor;
               else
                  cursor = predecessor(cursor);

            cursor = center;

            // recherche apres center
            while (cursor != null && ev.compareTo(cursor.event) == 0)
               if(cursor.event == ev)
                  return cursor;
               else
                  cursor = successor(cursor);

            return null;
         }
      }
      return null;
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

   /*
   public static void main (String[] args) {
      SplayTree sp = new SplayTree();
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
      sp.add(e5);
      sp.print2(sp.root);
      sp.add(e7);
      sp.print2(sp.root);
      sp.add(e5);
      sp.print2(sp.root);

      sp.getFirst();
      System.out.println(".....after Get" + PrintfFormat.NEWLINE +
                      PrintfFormat.NEWLINE + PrintfFormat.NEWLINE);
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
   }
*/
}
