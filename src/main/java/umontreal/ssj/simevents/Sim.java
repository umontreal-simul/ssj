/*
 * Class:        Sim
 * Description:  contains the executive of a discrete-event simulation
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
import umontreal.ssj.simevents.Event;
import umontreal.ssj.simevents.Simulator;

/**
 * This static class contains the executive of a discrete-event simulation.
 * It maintains the simulation clock and starts executing the events in the
 * appropriate order. Its methods permit one to start, stop, and
 * (re)initialize the simulation, and read the simulation clock.
 *
 * Starting from SSJ-2.0, the `Sim` class now uses the default simulator
 * returned by the `getDefaultSimulator()` method in the  @ref Simulator
 * class. Although the `Sim` class is perfectly adequate for simple
 * simulations, the  @ref Simulator class is more general and supports more
 * functionnalities. For example, if one needs to have more than one
 * simulation clock and event list, one will have to use the  @ref Simulator
 * class instead of the simpler `Sim` class.
 *
 * <div class="SSJ-bigskip"></div>
 */
public final class Sim {

   // Prevents construction of this object
   private Sim() {}

   /**
    * Returns the current value of the simulation clock.
    *  @return the current simulation time
    */
   public static double time() {
      return Simulator.getDefaultSimulator().time();
   }

   /**
    * Reinitializes the simulation executive by clearing up the event
    * list, and resetting the simulation clock to zero. This method must
    * not be used to initialize process-driven simulation;
    * umontreal.ssj.simprocs.SimProcess.init must be used instead.
    */
   public static void init() {
     // This has to be done another way in order to separate events and processes.
      Simulator.getDefaultSimulator().init();
   }

   /**
    * Same as  #init, but also chooses `evlist` as the event list to be
    * used. For example, calling `init(new DoublyLinked())` initializes
    * the simulation with a doubly linked linear structure for the event
    * list. This method must not be used to initialize process-driven
    * simulation;  umontreal.ssj.simprocs.DSOLProcessSimulator(init)
    * &nbsp;`(EventList)` or <br>
    * umontreal.ssj.simprocs.ThreadProcessSimulator(init)
    * &nbsp;`(EventList)` must be used instead.
    *  @param evlist       selected event list implementation
    */
   public static void init (EventList evlist) {
     Simulator.getDefaultSimulator().init(evlist);
   }

   /**
    * Gets the currently used event list.
    *  @return the currently used event list
    */
   public static EventList getEventList() {
      return Simulator.getDefaultSimulator().getEventList();
   }

   /**
    * This method is used by the package  @ref umontreal.ssj.simprocs; *it
    * should not be used directly by a simulation program*. It removes the
    * first event from the event list and sets the simulation clock to its
    * event time.
    *  @return the first planned event, or `null` if there is no such
    * event
    */
   protected static Event removeFirstEvent() {
       return Simulator.getDefaultSimulator().removeFirstEvent();
   }

   /**
    * Starts the simulation executive. There must be at least one `Event`
    * in the event list when this method is called.
    */
   public static void start() {
      Simulator.getDefaultSimulator().start();
   }

   /**
    * Tells the simulation executive to stop as soon as it takes control,
    * and to return control to the program that called  #start. This
    * program will then continue executing from the instructions right
    * after its call to `Sim.start`. If an  @ref Event is currently
    * executing (and this event has just called <tt>Sim.stop</tt>), the
    * executive will take control when the event terminates its execution.
    */
   public static void stop() {
      Simulator.getDefaultSimulator().stop();
   }


   // Used to passivate and activate the main process.
   // See the SimProcess.dispatch() and SimThread.actions()

//    protected static void activate() {
 //       synchronized (eventList) {eventList.notify(); }
 //    }

//     protected static void passivate() {
//         synchronized (eventList) {
  //          try { eventList.wait(); } catch (InterruptedException e) {}
  //       }
  //   }
}