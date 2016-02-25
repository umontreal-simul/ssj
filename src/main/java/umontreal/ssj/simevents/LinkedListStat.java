/*
 * Class:        LinkedListStat
 * Description:  
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
import java.util.LinkedList;
import java.util.NoSuchElementException;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class extends  @ref ListWithStat, and uses a linked list as the
 * internal data structure.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class LinkedListStat<E> extends ListWithStat<E> {

   /**
    * Constructs a new list, initially empty.
    */
   public LinkedListStat() {
      super (Simulator.getDefaultSimulator(), new LinkedList<Node<E>>());
   }

   /**
    * Constructs a new list, initially empty, and using the default
    * simulator.
    *  @param inSim        Simulator associate to the current variable.
    */
   public LinkedListStat(Simulator inSim) {
      super (inSim, new LinkedList<Node<E>>());
   }

   /**
    * Constructs a list containing the elements of the specified
    * collection, using the default simulator.
    *  @param c            collection containing elements to fill in this
    *                      list with
    */
   public LinkedListStat (Collection<? extends E> c) {
      super (Simulator.getDefaultSimulator(), new LinkedList<Node<E>>(), c);
   }

   /**
    * Constructs a list containing the elements of the specified
    * collection.
    *  @param inSim        Simulator associate to the current variable.
    *  @param c            collection containing elements to fill in this
    *                      list with
    */
   public LinkedListStat (Simulator inSim, Collection<? extends E> c) {
      super (inSim, new LinkedList<Node<E>>(), c);
   }

   /**
    * Constructs a new list with name `name`, using the default simulator.
    * This name can be used to identify the list in traces and reports.
    *  @param name         name for the list object
    */
   public LinkedListStat (String name) {
      super (Simulator.getDefaultSimulator(), new LinkedList<Node<E>>(), name);
   }

   /**
    * Constructs a new list with name `name`. This name can be used to
    * identify the list in traces and reports.
    *  @param inSim        Simulator associate to the current variable.
    *  @param name         name for the list object
    */
   public LinkedListStat (Simulator inSim, String name) {
      super (inSim, new LinkedList<Node<E>>(), name);
   }

   /**
    * Constructs a new list containing the elements of the specified
    * collection `c` and with name `name`, using the default simulator.
    * This name can be used to identify the list in traces and reports.
    *  @param c            collection containing elements to fill in this
    *                      list with
    *  @param name         name for the list object
    */
   public LinkedListStat (Collection<? extends E> c, String name) {
      super (Simulator.getDefaultSimulator(), new LinkedList<Node<E>>(), c, name);
   }

   /**
    * Constructs a new list containing the elements of the specified
    * collection `c` and with name `name`. This name can be used to
    * identify the list in traces and reports.
    *  @param inSim        Simulator associate to the current variable.
    *  @param c            collection containing elements to fill in this
    *                      list with
    *  @param name         name for the list object
    */
   public LinkedListStat (Simulator inSim, Collection<? extends E> c,
                          String name) {
      super (inSim, new LinkedList<Node<E>>(), c, name);
   }

   /**
    * @name `LinkedList` methods
    *
    * See the JDK documentation for more information about these methods.
    *
    * @{
    */
   public void addFirst (E obj) {
      add (0, obj);
   }
   public void addLast (E obj) {
      add (size(), obj);
   }
   public E getFirst() {
      if (isEmpty())
         throw new NoSuchElementException();
      return get (0);
    }
   public E getLast() {
      if (isEmpty())
         throw new NoSuchElementException();
      return get (size() - 1);
   }
   public E removeFirst() {
      if (isEmpty())
         throw new NoSuchElementException();
      return remove (0);
   }
   public E removeLast() {
      if (isEmpty())
         throw new NoSuchElementException();
      return remove (size() - 1);
   }
}

/**
 * @}
 */