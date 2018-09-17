/*
 * Class:        ListWithStat
 * Description:  Implements a list with integrated statistical probes to
                 provide automatic collection of statistics on the sojourn
                 times of objects in the list and on the size of the list
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
package umontreal.ssj.simevents;
import java.util.Collection;
import java.util.List;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.NoSuchElementException;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.TransformingList;
import umontreal.ssj.util.PrintfFormat;

/**
 * Implements a list with integrated statistical probes to provide automatic
 * collection of statistics on the sojourn times of objects in the list and
 * on the size of the list as a function of time given by a simulator. The
 * automatic statistical collection can be enabled or disabled for each list,
 * to reduce overhead. This class extends
 * @ref umontreal.ssj.util.TransformingList and transforms elements into
 * nodes associating insertion times with elements.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ListWithStat<E>
             extends TransformingList<E, ListWithStat.Node<E>> {
   private boolean stats; // true si on a appele setStatCollecting
   private double initTime; // temps de la derniere initialisation
   private Accumulate blockSize; //block stat. sur la longueur de la liste
   private Tally blockSojourn; // block stat. sur les durees de sejour
   private String name;
   private Simulator sim;

   /**
    * Constructs a new list with internal data structure using the default
    * simulator and implemented by `nodeList`. The given list is cleared
    * for the constructed list to be initially empty.
    */
   public ListWithStat (List<Node<E>> nodeList) {
      super (nodeList);
      nodeList.clear();
      sim = Simulator.getDefaultSimulator();
      stats = false;
   }

   /**
    * Constructs a new list with internal data structure implemented by
    * `nodeList`. The given list is cleared for the constructed list to be
    * initially empty.
    *  @param nodeList     the list containing the nodes
    */
   public ListWithStat (Simulator inSim, List<Node<E>> nodeList) {
      super (nodeList);
      if (inSim == null || nodeList == null)
          throw new NullPointerException();
      nodeList.clear();
      sim = inSim;
      stats = false;
   }

   /**
    * Constructs a list containing the elements of the specified
    * collection, whose elements are stored into `nodeList` and using the
    * default simulator.
    *  @param nodeList     the list containing the nodes
    *  @param c            collection containing elements to fill in this
    *                      list with
    */
   public ListWithStat (List<Node<E>> nodeList, Collection<? extends E> c) {
      this (Simulator.getDefaultSimulator(), nodeList);
      addAll (c);
   }

   /**
    * Constructs a list containing the elements of the specified
    * collection, whose elements are stored into `nodeList`.
    *  @param inSim        simulator associate to the current variable
    *  @param nodeList     the list containing the nodes
    *  @param c            collection containing elements to fill in this
    *                      list with
    */
   public ListWithStat (Simulator inSim, List<Node<E>> nodeList,
                        Collection<? extends E> c) {
      this (inSim, nodeList);
      addAll (c);
   }

   /**
    * Constructs a new list with name `name`, internal list `nodeList`,
    * and using the default simulator. This name can be used to identify
    * the list in traces and reports. The given list is cleared for the
    * constructed list to be initially empty.
    *  @param nodeList     the list containing the nodes
    *  @param name         name for the list object
    */
   public ListWithStat (List<Node<E>> nodeList, String name) {
      this (Simulator.getDefaultSimulator(), nodeList);
      this.name = name;
   }

   /**
    * Constructs a new list with name `name`, and internal list
    * `nodeList`. This name can be used to identify the list in traces and
    * reports. The given list is cleared for the constructed list to be
    * initially empty.
    *  @param inSim        simulator associate to the current variable
    *  @param nodeList     the list containing the nodes
    *  @param name         name for the list object
    */
   public ListWithStat (Simulator inSim, List<Node<E>> nodeList,
                        String name) {
      this (inSim, nodeList);
      this.name = name;
   }

   /**
    * Constructs a new list containing the elements of the specified
    * collection `c`, with name `name`, internal list `nodeList`, and
    * using the default simulator. This name can be used to identify the
    * list in traces and reports.
    *  @param nodeList     the list containing the nodes
    *  @param c            collection containing elements to fill in this
    *                      list with
    *  @param name         name for the list object
    */
   public ListWithStat (List<Node<E>> nodeList, Collection<? extends E> c,
                        String name) {
      this (Simulator.getDefaultSimulator(), nodeList);
      this.name = name;
      addAll (c);
   }

   /**
    * Constructs a new list containing the elements of the specified
    * collection `c`, with name `name`, and internal list `nodeList`. This
    * name can be used to identify the list in traces and reports.
    *  @param inSim        simulator associate to the current variable
    *  @param nodeList     the list containing the nodes
    *  @param c            collection containing elements to fill in this
    *                      list with
    *  @param name         name for the list object
    */
   public ListWithStat (Simulator inSim, List<Node<E>> nodeList,
                        Collection<? extends E> c, String name) {
      this (inSim, nodeList);
      this.name = name;
      addAll (c);
   }
public E convertFromInnerType (Node<E> node) {
      return node.getElement();
   }

   public Node<E> convertToInnerType (E element) {
      return new Node<E> (element, sim.time());
   }

/**
 * Returns the simulator associated with this list.
 *  @return the simulator associated with this list
 */
public Simulator simulator() {
      return sim;
   }

   /**
    * Sets the simulator associated with this list. This list should be
    * cleared after this method is called.
    *  @param sim          the simulator of this list
    */
   public void setSimulator(Simulator sim) {
       if (sim == null)
         throw new NullPointerException();
      this.sim = sim;
      if (blockSize != null)
         blockSize.setSimulator (sim);
   }


   @Override
   public void clear() {
      if (stats)
         initStat();
      super.clear();
   }

   @Override
   public void add (int index, E obj)  {
      super.add (index, obj);
      if (stats)
         blockSize.update (size());
   }

   @Override
   public E remove (int index) {
      Node<E> node = getInnerList().get (index);
      if (stats)
         blockSojourn.add (sim.time() - node.getInsertionTime());
      E e = super.remove (index);
      if (stats)
         blockSize.update (size());
      return e;
   }

   @Override
   public Iterator<E> iterator() {
      return new IteratorWithStat (getInnerList().iterator());
   }

   @Override
   public ListIterator<E> listIterator() {
      return new ListIteratorWithStat (getInnerList().listIterator());
   }

   @Override
   public ListIterator<E> listIterator (int index) {
      return new ListIteratorWithStat (getInnerList().listIterator (index));
   }

   @Override
   public E set (int index, E element) {
      Node<E> oldNode = getInnerList().get (index);
      E oldElement = oldNode.getElement();
      boolean equal;
      if (oldElement == null || element == null)
         equal = oldElement == element;
      else
         equal = oldElement.equals (element);
      if (equal) {
         getInnerList().set (index, new Node<E> (element, oldNode.getInsertionTime()));
         return oldElement;
      }
      else {
         if (stats)
            blockSojourn.add (sim.time() - oldNode.getInsertionTime());
         getInnerList().set (index, new Node<E> (element, sim.time()));
         return oldElement;
      }
   }

   private class IteratorWithStat implements Iterator<E> {
      private Iterator<Node<E>> itr;
      private Node<E> lastRet;

      public IteratorWithStat (Iterator<Node<E>> itr) {
         this.itr = itr;
      }

      public boolean hasNext() {
         return itr.hasNext();
      }

      public E next() {
         lastRet = itr.next();
         return lastRet.getElement();
      }

      public void remove() {
         itr.remove();
         if (stats) {
            blockSize.update (size());
            blockSojourn.add (sim.time() - lastRet.getInsertionTime());
         }
         lastRet = null;
      }
   }

   private class ListIteratorWithStat implements ListIterator<E> {
      private ListIterator<Node<E>> itr;
      private Node<E> lastRet;

      public ListIteratorWithStat (ListIterator<Node<E>> itr) {
         this.itr = itr;
      }

      public void add (E o) {
         itr.add (new Node<E> (o, sim.time()));
         lastRet = null;
         if (stats)
            blockSize.update (size());
      }

      public boolean hasNext() {
         return itr.hasNext();
      }

      public boolean hasPrevious() {
         return itr.hasPrevious();
      }

      public E next() {
         lastRet = itr.next();
         return lastRet.getElement();
      }

      public int nextIndex() {
         return itr.nextIndex();
      }

      public E previous() {
         lastRet = itr.previous();
         return lastRet.getElement();
      }

      public int previousIndex() {
         return itr.previousIndex();
      }

      public void remove() {
         itr.remove();
         if (stats) {
            blockSize.update (size());
            blockSojourn.add (sim.time() - lastRet.getInsertionTime());
         }
         lastRet = null;
      }

      public void set (E element) {
         if (lastRet == null)
            throw new NoSuchElementException();
         Node<E> oldNode = lastRet;
         E oldElement = oldNode.getElement();
         boolean equal;
         if (oldElement == null || element == null)
            equal = oldElement == element;
         else
            equal = oldElement.equals (element);
         if (equal) {
            lastRet = new Node<E> (element, oldNode.getInsertionTime());
            itr.set (lastRet);
         }
         else {
            if (stats)
               blockSojourn.add (sim.time() - oldNode.getInsertionTime());
            lastRet = new Node<E> (element, sim.time());
            itr.set (lastRet);
         }
      }
   }

  /**
   * @name Statistic collection methods
   * @{
   */

  /**
   * Returns `true` if the list collects statistics about its size and
   * sojourn times of elements, and `false` otherwise. By default,
   * statistical collecting is turned off.
   *  @return the status of statistical collecting
   */
  public boolean getStatCollecting() {
     return stats;
  }

  /**
   * Starts or stops collecting statistics on this list. If the statistical
   * collection is turned ON, the method creates two statistical probes if
   * they do not exist yet. The first one, of the class  @ref Accumulate,
   * measures the evolution of the size of the list as a function of time.
   * It can be accessed by the method  #statSize. The second one, of the
   * class  @ref umontreal.ssj.stat.Tally and accessible via  #statSojourn,
   * samples the sojourn times in the list of the objects removed during
   * the observation period, i.e., between the last initialization time of
   * this statistical probe and the current time. The method automatically
   * calls  #initStat to initialize these two probes. When this method is
   * used, it is normally invoked immediately after calling the constructor
   * of the list.
   *  @exception IllegalStateException if the statistical collection is in
   * the same state as the caller requires
   */
  public void setStatCollecting (boolean b) {
    if (b && !stats) {
      if (blockSize == null)
      blockSize = new Accumulate(sim, "List Size " + name);
      if (blockSojourn == null)
      blockSojourn = new Tally("List Sojourn " + name);
      blockSize.update (size());
      stats = true;
      initStat();
    } else
    stats = false;
  }

   /**
    * Reinitializes the two statistical probes created by
    * {@link #setStatCollecting() setStatCollecting(true)} and makes an
    * update for the probe on the list size.
    *  @exception IllegalStateException if the statistical collection is
    * disabled
    */
   public void initStat() {
    if (!stats)
    throw new IllegalStateException("initStat for a list that did not call setStatCollecting (true).");
    blockSize.init();
    blockSojourn.init();
    blockSize.update (size());
    initTime = sim.time();
  }

   /**
    * Returns the last simulation time  #initStat was called.
    *  @return the last simulation time  #initStat was called
    */
   public double getInitTime() {
      return initTime;
   }

   /**
    * Returns the statistical probe on the evolution of the size of the
    * list as a function of the simulation time. This probe exists only if
    * {@link #setStatCollecting() setStatCollecting(true)} has been called
    * for this list.
    *  @return the statistical probe on the evolution of the size of the
    * list
    */
   public Accumulate statSize()   {
       return blockSize;
   }

   /**
    * Returns the statistical probe on the sojourn times of the objects in
    * the list. This probe exists only if  {@link #setStatCollecting()
    * setStatCollecting(true)} has been called for this list.
    *  @return the statistical probe for the sojourn times in the list
    */
   public Tally statSojourn()  {
      return blockSojourn;
   }

   /**
    * Returns a string containing a statistical report on the list,
    * provided that  {@link #setStatCollecting() setStatCollecting(true)}
    * has been called before for this list. Even If  #setStatCollecting
    * was called with `false` afterward, the report will be made for the
    * collected observations. If the probes do not exist, i.e.,
    * #setStatCollecting was never called for this object, an illegal
    * state exception will be thrown.
    *  @return a statistical report, represented as a string
    *
    *  @exception IllegalStateException if no statistical probes exist
    */
   public String report()  {
        if (blockSojourn == null || blockSize == null)
            throw new IllegalStateException
                ("Calling report when no statistics were collected");

        PrintfFormat str = new PrintfFormat();
        str.append (PrintfFormat.NEWLINE +
            "REPORT ON LIST : ").append (name).append (PrintfFormat.NEWLINE);
        str.append ("   From time: ").append (7, 2, 2, initTime);
        str.append (" to time: ").append (10, 2, 2, sim.time());
        str.append ("                  min        max      average  ");
        str.append ("standard dev.  nb. Obs");

        str.append ("   Size    ");
        str.append (9, (int)(blockSize.min()+0.5));
        str.append (11, (int)(blockSize.max()+0.5));
        str.append (14, 3, 2, blockSize.average()).append(PrintfFormat.NEWLINE);

        str.append ("   Sojourn ");
        str.append ( 12, 3, 2, blockSojourn.min()).append (" ");
        str.append (10, 3, 2, blockSojourn.max()).append (" ");
        str.append (10, 3, 2, blockSojourn.average()).append (" ");
        str.append (10, 3, 2, blockSojourn.standardDeviation()).append (" ");
        str.append (11, blockSojourn.numberObs()).append (PrintfFormat.NEWLINE);

        return str.toString();
    }

   /**
    * Returns the name associated to this list, or `null` if no name was
    * assigned.
    *  @return the name associated to this list
    */
   public String getName() {
      return name;
   }

   /**
    * @}
    */

   /**
    * @name Inner class
    * @{
    */

   /**
    * Represents a node that can be part of a list with statistical
    * collecting.
    */
   public static class Node<E> {
      private E element;
      private double insertionTime;

      /**
       * @}
       */

      /**
       * Constructs a new node containing element `element` inserted
       * into the list at time `insertionTime`.
       *  @param element      the element to add into this new node
       *  @param insertionTime the insertion time of the element
       */
      public Node (E element, double insertionTime) {
         this.element = element;
         this.insertionTime = insertionTime;
      }

      /**
       * Returns the element stored into this node.
       *  @return the element into this node
       */
      public E getElement() { return element; }

      /**
       * Returns the insertion time of the element in this node.
       *  @return the insertion time of the element
       */
      public double getInsertionTime() { return insertionTime; }


      public String toString() {
         String str = element == null ? "null" : element.toString();
         str += " (inserted at time " + insertionTime + ")";
         return str;
      }
   }

}