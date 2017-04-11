/*
 * Class:        Henriksen
 * Description:  implementation of class EventList using the doubly-linked
                 indexed list of Henriksen
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
 * An implementation of  @ref EventList using the doubly-linked indexed list
 * of Henriksen @cite sKIN86a&thinsp; (see also @cite sFIS01a&thinsp; (p.
 * 207)).
 *
 * Events are stored in a normal doubly-linked list. An additional index
 * array is added to the structure to allow quick access to the events.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Henriksen implements EventList {
   /*
     Fonctionnement de l'algorithme :

     On maintient une liste doublement chainee contenant les evenements.
     Les deux extremites de cette liste sont occupees par des bornes qui
     ont des temps qui ne peuvent pas etre depasses.

     Au-dessus de cette liste se retrouve un index (trie en ordre decroissant)
     qui contient des pointeurs vers certaines des entrees de la liste
     chainee. Le premier element de l'index est toujours la borne superieure.
     On se sert de cet index pour faire des recherches binaires pour
     retrouver rapidement un element dans la liste. On fait d'abord une
     recherche binaire dans l'index avant de faire une recherche lineaire
     dans la liste. On compte aussi le nombre d'elements qui doivent etre
     parcourus dans la recherche lineaire. Si ce nombre atteint une certaine
     valeur (4 dans cette implantation), alors on fait en sorte qu'un
     element de l'index pointe vers l'element de la liste que l'on parcourt.
     Ceci permet un certain balancement des elements presents dans l'index,
     et donc une meilleure efficacite de la recherche binaire sans avoir
     a mettre tous les elements dans l'index.
   */

   private static final double MIN_VALUE = -10E38; //n'importe quoi < 0
   private static final double MAX_VALUE = 10E38;  //le plus gros possible

   private static final int ARRAY_LENGTH_INIT = 256;


   private int modCount = 0;

   private Entry firstEntry;

   //for the binary search
   private Entry[] entryVec;

   private int vectSize;
   private int arrayLength;


   public Henriksen() {
      //creation des bornes
      Entry lastEntry = new Entry(null, null, null, MAX_VALUE);
      firstEntry = new Entry(null, null, lastEntry, MIN_VALUE);
      lastEntry.left = firstEntry;

      arrayLength = ARRAY_LENGTH_INIT;
      entryVec = new Entry[arrayLength];

      changeSize(1);
      entryVec[0] = lastEntry;
   }

   @Override
   public boolean isEmpty() {
      return firstEntry.right == entryVec[0];
   }

   @Override
   public void clear() {
      if(isEmpty())
         return;

      firstEntry.right = entryVec[0];
      entryVec[0].left = firstEntry;
      changeSize(1);

      //on enleve tous les liens menant aux entrees supprimees :
      for(int i = 1; i < arrayLength; i++)
         entryVec[i] = null;

      modCount++;
   }


   @Override
   public void add (Event ev) {
      Entry prec = findEntry(ev, false);

      Entry e = new Entry(ev, prec, prec.right, ev.time());
      e.right.left = e;
      prec.right = e;
      modCount++;
   }

   @Override
   public void addFirst (Event ev) {
      Entry e = new Entry(ev, firstEntry, firstEntry.right, ev.time());
      firstEntry.right.left = e;
      firstEntry.right = e;

      modCount++;
   }

   @Override
   public void addBefore (Event ev, Event other) {
      Entry otherEntry = findEntry(other, true);
      if(otherEntry == null)
         throw new IllegalArgumentException("Event not in list.");
      Entry e = new Entry(ev, otherEntry.left, otherEntry, ev.time());
      otherEntry.left.right = e;
      otherEntry.left = e;

      modCount++;
   }

   @Override
   public void addAfter (Event ev, Event other) {
      Entry otherEntry = findEntry(other, true);
      if(otherEntry == null)
         throw new IllegalArgumentException("Event not in list.");
      Entry e = new Entry(ev, otherEntry, otherEntry.right, ev.time());
      otherEntry.right.left = e;
      otherEntry.right = e;

      modCount++;
   }


   @Override
   public Event getFirst() {
      return firstEntry.right.event;
   }

   @Override
   public Event getFirstOfClass (String cl) {
      Entry e = firstEntry.right;
      while(e.right != null) {
         if(e.event.getClass().getName().equals(cl))
            return e.event;
         e = e.right;
      }
      return null;
   }

   @SuppressWarnings("unchecked")
   @Override
   public <E extends Event> E getFirstOfClass (Class<E> cl) {
      Entry e = firstEntry.right;
      while(e.right != null) {
         if(e.event.getClass() == cl)
            return (E)e.event;
         e = e.right;
      }
      return null;
   }

   @Override
   public Iterator<Event> iterator() {
      return listIterator();
   }

   @Override
   public ListIterator<Event> listIterator() {
      return new HItr();
   }

   @Override
   public boolean remove (Event ev) {
      Entry e = findEntry (ev, true);
      if (e == null)
         return false;

      //on l'enleve de l'index
      int i = findIndex (ev.time());
      i++;

      while (i < vectSize && entryVec[i].event != null &&
             ev.compareTo(entryVec[i].event) == 0) {
         if(entryVec[i].event == ev)
            entryVec[i] = e.left;

         i++;
      }

      //on l'enleve de la liste
      e.left.right = e.right;
      e.right.left = e.left;
      e.right = null;
      e.left = null;
      e.event = null;

      modCount++;

      return true;
   }

   @Override
   public Event removeFirst() {
      // si la premiere moitie de l'index est composee d'entrees perimees,
      // on coupe de moitie l'index
      if (entryVec[vectSize/2].time <= firstEntry.right.time && vectSize > 1)
         changeSize(vectSize / 2);

      Entry e = firstEntry.right;

      //borne superieure
      if (e == entryVec[0])
         return null;

      firstEntry.right = e.right;
      e.right.left = firstEntry;

      e.right = null;
      e.left = null;

      Event ev = e.event;
      e.event = null;

      modCount++;

      return ev;
   }


   @Override
   public String toString() {
      StringBuilder sb = new StringBuilder
                        ("Contents of the event list Henriksen:");
      Entry e = firstEntry.right;
      while (e.right != null) {
         sb.append (PrintfFormat.NEWLINE +
                    PrintfFormat.g (12, 7, e.event.time()) + ", " +
                    PrintfFormat.g (8, 4, e.event.priority()) + " : " +
                    e.event.toString());
         e = e.right;
      }
      return sb.toString();
   }


   private static class Entry {
      public Event event;
      public Entry left;
      public Entry right;
      public double time;

      Entry (Event event, Entry left, Entry right, double time) {
         this.event = event;
         this.left = left;
         this.right = right;
         this.time = time;
      }

      @Override
      public String toString() {
         return "[" + event + " |" + time + "|]";
      }
   }


   private class HItr implements ListIterator<Event> {
      private Entry prev;
      private Entry next;
      private Entry lastRet;
      private int expectedModCount;
      private int nextIndex;

      private HItr() {
         prev = firstEntry;
         next = firstEntry.right;
         expectedModCount = modCount;
         lastRet = null;
         nextIndex = 0;
      }

      @Override
      public void add (Event ev) {
         if(modCount != expectedModCount)
            throw new ConcurrentModificationException();

         //make sure the time is in the right order
         if (ev.time() > next.time) {
            ev.setTime (next.time);
            ev.setPriority (next.event.priority());
         }
         if (ev.time() < prev.time) {
            ev.setTime (prev.time);
            ev.setPriority (prev.event.priority());
         }

         Entry e = new Entry(ev, prev, next, ev.time());
         prev.right = e;
         next.left = e;
         prev = e;

         nextIndex++;
         lastRet = null;
         modCount++;
         expectedModCount++;
      }

      @Override
      public boolean hasNext() {
         if(modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return next != entryVec[0];
      }

      @Override
      public boolean hasPrevious() {
         if(modCount != expectedModCount)
            throw new ConcurrentModificationException();
         return next != firstEntry;
      }

      @Override
      public Event next() {
         if (!hasNext())
            throw new NoSuchElementException();

         nextIndex++;
         Event ev = next.event;
         lastRet = next;
         next = next.right;
         prev = prev.right;
         return ev;
      }

      @Override
      public int nextIndex() {
         if(!hasNext())
            throw new NoSuchElementException();

         return nextIndex;
      }

      @Override
      public Event previous() {
         if(!hasPrevious())
            throw new NoSuchElementException();

         nextIndex--;
         Event ev = prev.event;
         lastRet = prev;
         prev = prev.left;
         next = next.left;
         return ev;
      }

      @Override
      public int previousIndex() {
         if(!hasPrevious())
            throw new NoSuchElementException();

         return nextIndex - 1;
      }


      @Override
      public void remove() {
         if(modCount != expectedModCount)
            throw new ConcurrentModificationException();
         if(lastRet == null)
            throw new IllegalStateException();

         if(lastRet == next) { //last call was to previous, not next
            if(next != entryVec[0])
               next = next.right;
         } else { //last call was to next or nothing
            if(prev != firstEntry) {
               prev = prev.left;
               nextIndex--;
            }
         }

         //remove the deleted entry in the vector
         double evtime = lastRet.time;
         int i = findIndex (evtime);
         i++;

         while(i < vectSize && entryVec[i].time == evtime) {
            if(entryVec[i].event == lastRet.event)
               entryVec[i] = lastRet.left;
            i++;
         }

         lastRet.event = null;
         lastRet.left.right = lastRet.right;
         lastRet.right.left = lastRet.left;
         lastRet.left = null;
         lastRet.right = null;
         lastRet = null;
         modCount++;
         expectedModCount++;
      }

      @Override
      public void set (Event ev) {
         if(modCount != expectedModCount)
            throw new ConcurrentModificationException();
         if(lastRet == null)
            throw new IllegalStateException();

         // Check for a good time
         if (ev.time() < lastRet.left.time) {
            ev.setTime (lastRet.left.time);
            ev.setPriority (lastRet.left.event.priority());
         }
         if (ev.time() > lastRet.right.time) {
            ev.setTime (lastRet.right.time);
            ev.setPriority (lastRet.right.event.priority());
         }

         lastRet.event = ev;
      }
   }

   /*
     On change la taille de l'index.
     Le changement de taille se fait du cote des entrees de temps inferieures.
     Ces entrees se retrouvent a la fin de l'index pour simplifier ce
     travail.
    */
   private void changeSize(int newSize) {
      // si on grossit le vecteur reel
      if(newSize > arrayLength) {
         Entry[] newVec = new Entry[newSize];
         for(int i = 0; i < vectSize; i++)
            newVec[i] = entryVec[i];
         entryVec = newVec;
         arrayLength = newSize;
      }

      // les nouveaux emplacements sont remplis par la borne min
      for(int i = vectSize; i < newSize; i++)
         entryVec[i] = firstEntry;

      vectSize = newSize;
   }

   /*
     Fait une recherche binaire a l'interieur de l'index pour trouve
     l'entree dans l'index qui a la plus petite valeur de temps, mais
     dont la valeur est superieure a evtime.
     Les entrees sont triees en ordre inverse dans l'index, pour simplifie
     le changement de taille de l'index.
    */
   private int findIndex (double evtime) {
      int i = vectSize / 2;
      int j = vectSize / 4;

      // recherche binaire dans le vecteur d'index
      while(j > 0) {
         // note : entryVec est trie a l'envers
         if(evtime >= entryVec[i].time)
            i -= j;
         else
            i += j;
         j /= 2;
      }

      if (evtime >= entryVec[i].time)
         i--;

      return i;
   }

   /*
     Si findEvent est false, trouve la derniere entree qui a un temps
     egal ou inferieur au temps de l'evenement ev.
     Sinon, trouve l'entree contenant ev dans la liste.
    */
   private Entry findEntry (Event ev, boolean findEvent) {
      double evtime = ev.time();
      int i = findIndex(evtime);

      Entry e = entryVec[i].left;
      if (null == e) return null;
      int count = 0;

      while (e.time >= evtime/* && e.event != null*/ && ev.compareTo(e.event) < 0) {
         ++count;
         if (count == 4) {
            //operation pull
            if (i+1 >= vectSize)
               changeSize(vectSize * 2);

            i++;
            count = 0;
            entryVec[i] = e;
         }

         e = e.left;
      }

      if (findEvent) {
         // on cherche l'evenement identique
         Entry start = e;

         while (e != firstEntry && e.time == evtime && e.event != ev)
            e = e.left;

         // on ne l'a pas trouve
         if (e.event != ev)
            return null;
      }

      return e;
   }

}