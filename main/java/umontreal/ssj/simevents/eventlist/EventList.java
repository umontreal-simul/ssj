/*
 * Class:        EventList
 * Description:  interface for implementations of event lists
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

import java.util.ListIterator; 
import umontreal.ssj.simevents.Event;

/**
 * An interface for implementations of event lists. Different implementations
 * are provided in SSJ: doubly-linked list, splay tree, Henricksen’s method,
 * etc. The *events* in the event list are objects of the class
 * @ref umontreal.ssj.simevents.Event. The method
 * umontreal.ssj.simevents.Sim.init(EventList) permits one to select the
 * actual implementation used in a simulation @cite sKIN85a&thinsp;.
 *
 * To allow the user to print the event list, the  java.lang.Object.toString
 * method from the  Object class should be reimplemented in all
 * <tt>EventList</tt> implementations. It will return a string in the
 * following format: "<tt>Contents of the event list </tt><em>event list
 * class</em><tt>:</tt>" for the first line and each subsequent line has
 * format "<em>scheduled event time</em><tt>, </tt><em>event
 * priority</em>&nbsp;<tt>:</tt>&nbsp;<em>event string</em>". The *event
 * string* is obtained by calling the <tt>toString</tt> method of the event
 * objects. The string should not end with the end-of-line character.
 *
 * The following example is the event list of the bank example, printed at
 * 10h30. See <tt>examples.pdf</tt> for more information.
 *
 * <tt>
 * <pre>
 * Contents of the event list SplayTree:
 *    10.51,        1 : BankEv$Arrival\@cfb549
 *    10.54,        1 : BankEv$Departure\@8a7efd
 *       11,        1 : BankEv$3\@86d4c1
 *       14,        1 : BankEv$4\@f9f9d8
 *       15,        1 : BankEv$5\@820dda
 * </pre>
 * </tt>
 *
 * <div class="SSJ-bigskip"></div>
 */
public interface EventList extends Iterable<Event> {

/**
 * Returns <tt>true</tt> if and only if the event list is empty (no event is
 * scheduled).
 *  @return <tt>true</tt> if the event list is empty
 */
public boolean isEmpty();

   /**
    * Empties the event list, i.e., cancels all events.
    */
   public void clear();

   /**
    * Adds a new event in the event list, according to the time of
    * <tt>ev</tt>. If the event list contains events scheduled to happen
    * at the same time as <tt>ev</tt>, <tt>ev</tt> must be added after all
    * these events.
    *  @param ev           event to be added
    */
   public void add (Event ev);

   /**
    * Adds a new event at the beginning of the event list. The given event
    * <tt>ev</tt> will occur at the current simulation time.
    *  @param ev           event to be added
    */
   public void addFirst (Event ev);

   /**
    * Same as  #add, but adds the new event <tt>ev</tt> immediately before
    * the event <tt>other</tt> in the list.
    *  @param ev           event to be added
    *  @param other        reference event before which <tt>ev</tt> will
    *                      be added
    */
   public void addBefore (Event ev, Event other);

   /**
    * Same as  #add, but adds the new event <tt>ev</tt> immediately after
    * the event <tt>other</tt> in the list.
    *  @param ev           event to be added
    *  @param other        reference event after which <tt>ev</tt> will be
    *                      added
    */
   public void addAfter (Event ev, Event other);

   /**
    * Returns the first event in the event list. If the event list is
    * empty, returns <tt>null</tt>.
    *  @return the first event in the event list, or <tt>null</tt> if the
    * list is empty
    */
   public Event getFirst();

   /**
    * Returns the first event of the class <tt>cl</tt> (a subclass of
    * <tt>Event</tt>) in the event list. If no such event is found,
    * returns <tt>null</tt>.
    *  @return the first event of class <tt>cl</tt>, or <tt>null</tt> if
    * no such event exists in the list
    */
   public Event getFirstOfClass (String cl);

   /**
    * Returns the first event of the class <tt>E</tt> (a subclass of
    * <tt>Event</tt>) in the event list. If no such event is found,
    * returns <tt>null</tt>.
    *  @return the first event of class <tt>cl</tt>, or <tt>null</tt> if
    * no such event exists in the list
    */
   public <E extends Event> E getFirstOfClass (Class<E> cl);

   /**
    * Returns a list iterator over the elements of the class
    * <tt>Event</tt> in this list.
    *  @return a list iterator over the elements of the class
    * <tt>Event</tt> in this list
    */
   public ListIterator<Event> listIterator();

   /**
    * Removes the event <tt>ev</tt> from the event list (cancels this
    * event). Returns <tt>true</tt> if and only if the event removal has
    * succeeded.
    *  @param ev           event to be removed
    *  @return <tt>true</tt> if the event was successfully removed from
    * the list
    */
   public boolean remove (Event ev);

   /**
    * Removes the first event from the event list (to cancel or execute
    * this event). Returns the removed event. If the list is empty, then
    * <tt>null</tt> is returned.
    *  @return the first event removed from the list, or <tt>null</tt> if
    * the list is empty
    */
   public Event removeFirst();

}