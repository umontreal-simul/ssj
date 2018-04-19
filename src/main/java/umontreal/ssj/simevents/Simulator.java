/*
 * Class:        Simulator
 * Description:  Represents the executive of a discrete-event simulator
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

import umontreal.ssj.simevents.eventlist.EventList;
import umontreal.ssj.simevents.eventlist.SplayTree;

/**
 * Represents the executive of a discrete-event simulator. This class
 * maintains a simulation clock, an event list, and starts executing the
 * events in the appropriate order. Its methods permit one to start, stop,
 * and (re)initialize the simulation, and read the simulation clock.
 *
 * Usually, a simulation program uses a single simulation clock which is
 * represented by an instance of this class. For more convenience and
 * compatibility, this class therefore provides a mechanism to construct and
 * return a default simulator which is used when an event is constructed
 * without an explicit reference to a simulator, and when the simulator is
 * accessed through the  @ref Sim class.
 *
 * Note that this class is NOT thread-safe. Consequently, if a simulation
 * program uses multiple threads, it should acquire a lock on a simulator
 * (using a `synchronized` block) before accessing its state. Note however,
 * that one can launch many simulations in parallel with as many threads, as
 * long as *each thread has its own* `Simulator`.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Simulator {

   protected double currentTime = 0.0;
      // The current simulation time (clock).

   protected EventList eventList;
      // The list of future events.
      // Can be changed by the method \texttt{init}.

   protected boolean stopped = true;
      // Becomes true when the simulation has ended
      // (stopped has been called or the event list is empty).

   protected boolean simulating = false;

   protected ContinuousState continuousState = null;

   /**
    * Represents the default simulator being used by the class  @ref Sim, and
    * the no-argument constructor of  @ref Event. This simulator is usually
    * obtained with the  #getDefaultSimulator method, which initializes it if
    * needed. But it might also be initialized differently, e.g., if
    * process-driven simulation is required.
    */
   public static Simulator defaultSimulator = null;

   /**
    * Constructs a new simulator using a splay tree for the event list.
    */
   public Simulator() {
     eventList  = new SplayTree();
   }

   /**
    * Constructs a new simulator using `eventList` for the event list.
    */
   public Simulator (EventList eventList) {
     if (eventList == null)
        throw new NullPointerException();
     this.eventList = eventList;
   }

   /**
    * Returns the current value of the simulation clock.
    *  @return the current simulation time
    */
   public double time() {
      return currentTime;
   }

   /**
    * Reinitializes the simulation executive by clearing up the event
    * list, and resetting the simulation clock to zero.
    */
   public void init() {
     // This has to be done another way in order to separate events and processes.
//      SimProcess.killAll();
      currentTime = 0.0;   eventList.clear();   stopped = false;  simulating = false;
   }

   /**
    * Same as  #init, but also sets `evlist` as the event list to be
    * used. For example, `init(new DoublyLinked())` initializes the
    * simulator with a doubly linked linear structure 
    * @ref umontreal.ssj.simevents.eventlist.DoublyLinked for the event list.
    * Note that this method will clear the event list, so any 
    * existing events in `evlist` will be removed.
    *
    * To initialize the simulator with a non-empty event list,
    * one can extend the class `Simulator` and define his own initialization method.
    *
    *  @param evlist       selected event list implementation
    */
   public void init (EventList evlist) {
      if (evlist == null)
         throw new NullPointerException();
      eventList = evlist;
      init();  // will clear the events in evlist
   }

   /**
    * Gets the currently used event list.
    *  @return the currently used event list
    */
   public EventList getEventList() {
      return eventList;
   }

   /**
    * Determines if this simulator is currently running, i.e., executing
    * scheduled events.
    */
   public boolean isSimulating() {
       return simulating;
   }

   /**
    * Determines if this simulator was stopped by an event. The simulator
    * may still be processing the event which has called the  #stop
    * method; in this case,  #isSimulating returns `true`.
    */
   public boolean isStopped() {
       return stopped;
   }

   /**
    * Removes the first event from the event list and sets the simulation
    * clock to its event time.
    *  @return the first planned event, or `null` if there is no such
    * event
    */
   protected Event removeFirstEvent() {
       if (stopped)
          return null;
       Event ev = eventList.removeFirst();
       if (ev == null)
          return null;
       currentTime = ev.eventTime;
       ev.eventTime = -10.0;
       return ev;
   }

   /**
    * Starts the simulation executive. There must be at least one `Event`
    * in the event list when this method is called.
    */
   public void start () {
      if (eventList.isEmpty())
        throw new IllegalStateException ("start() called with an empty event list");
      stopped = false;
      simulating = true;
      Event ev;
      try {
         while ((ev = removeFirstEvent()) != null && !stopped) {
   //      while (!stopped && (ev = eventList.removeFirst()) != null) {
   //          currentTime = ev.eventTime;
   //          ev.eventTime = -10.0;
             ev.actions();
             // if ev is a thread object associated to a process,
             // the control will be transfered to this thread and the
             // executive will be passivated in the actions() method.
         }
      }
      finally {
         stopped = true; simulating = false;
      }
   }

   /**
    * Tells the simulation executive to stop as soon as it takes control,
    * and to return control to the program that called  #start. This
    * program will then continue executing from the instructions right
    * after its call to  #start(.) If an  @ref Event is currently
    * executing (and this event has just called  #stop ), the executive
    * will take control when the event terminates its execution.
    */
   public void stop()
   {
      stopped = true;
   }

   /**
    * Returns the current state of continuous variables being integrated
    * during the simulation. This state is used by the  @ref Continuous
    * class when performing simulation of continuous variables; it
    * defaults to an empty state, which is initialized only when this
    * method is called.
    *  @return continuousState field
    */
   public ContinuousState continuousState()
   {
      if (continuousState == null)
         continuousState = new ContinuousState(this);
      return continuousState;
   }

   /**
    * @name Static methods
    * @{
    */

   /**
    * Returns the default simulator instance used by the deprecated class
    * @ref Sim. If this simulator does not exist yet, it is constructed
    * using the no-argument constructor of this class. One can specify a
    * different default simulator by setting the `defaultSimulator` field
    * directly.
    */
   public static Simulator getDefaultSimulator() {
      if (defaultSimulator == null)
            defaultSimulator = new Simulator();
      return defaultSimulator;
   }

}

/**
 * @}
 */
