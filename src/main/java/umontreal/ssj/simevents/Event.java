/*
 * Class:        Event
 * Description:  provides event scheduling tools
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

/**
 * This abstract class provides event scheduling tools. Each type of event
 * should be defined as a subclass of the class `Event`, and should provide
 * an implementation of the method  #actions which is executed when an event
 * of this type occurs. The instances of these subclasses are the actual
 * events.
 *
 * Each event is linked to a simulator represented by an instance of
 * @ref Simulator before it can be scheduled and processed. A default
 * simulator, given by `Simulator.getDefaultSimulator`, is used if no
 * simulator is linked explicitly with an event. When an event is
 * constructed, it is not scheduled. It must be scheduled separately by
 * calling one of the scheduling methods  #schedule,  #scheduleNext,
 * #scheduleBefore, etc. An event can also be cancelled before it occurs.
 *
 * A scheduled event has an associated time at which it will happen and a
 * priority, which can be queried using the methods  #time and  #priority,
 * respectively. By default, events occur in ascending order of time, and
 * have priority 1. Events with the same time occur in ascending order of
 * priority. For example, if events `e1` and `e2` occur at the same time with
 * priority 2 and 1 respectively, then `e2` will occur before `e1`. Events
 * with the same time and priority occur in the order they were scheduled.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class Event implements Comparable<Event> {

   protected Simulator sim;
   //simulator linked with the current event

   protected double priority;
   //priority of the event. Priority is a second parameter (after eventTime)
   // used to class events for their running order, in the EventList.

   protected double eventTime;
   // Planned time of occurence of this event.  Negative if not planned.
   // Is protected because it is used (changed) in Process.

   // Replace that with instanceof simProcess to completely detach processes.
   // protected boolean isProcess = false;
   // Will be true for objects of the subclass Process of the class Event.
   // (i.e., true if this event is a process.)

   private int myra = 0;
   // A new event must always occur after those with same time and 
   // same priority in the Event list. myra is used for that in
   // SplayTree.java.

   // For internal use
   public final int getRa() { return myra; }
   public final void setRa(int r) { myra = r; }

   /**
    * Constructs a new event instance, which can be placed afterwards into
    * the event list of the default simulator by calling one of the
    * `schedule...` variants. For example, if `Bang` is an `Event`
    * subclass, the statement "<tt>new Bang().scheduleNext();</tt>"
    * creates a new `Bang` event and places it at the beginning of the
    * event list.
    * @code
    *
    *    public Event (double delay) \begin{hide} {
    *       this (Simulator.getDefaultSimulator(), delay);
    *    } \end{hide}
    *
    *    public Event (Simulator sim, double delay)\begin{hide} {
    *       if (sim == null)
    *          throw new NullPointerException();
    *       if (delay >= 0.0) {
    *          priority = 1;
    *          this.sim = sim;
    *          eventTime = sim.time() + delay;
    *          sim.eventList.add (this);
    *       }
    *       else
    *          throw new IllegalArgumentException ("Cannot schedule in the
    * past.");
    *    }\end{hide}
    *
    * @endcode
    *  Constructs a new event and inserts it in the event list of the
    * default simulator. If `delay >= 0.0`, the event is scheduled to
    * happen in `delay` units of simutated time. If two or more events are
    * scheduled to happen at the same time, events with the highest
    * priorities (lowest value of the `priority` field) occur first. If
    * two or more events are schedule to the same time, with the same
    * priority, they are placed in the event list (and executed) in the
    * same order as they have been scheduled.
    *
    * We recall that such constructors with parameters are not inherited
    * automatically by the subclasses in Java, but they can be invoked
    * using `super`. For example, one can have
    *
    * @code
    *
    *    class Bang extends Event {
    *       public Bang (double delay) { super (delay); }
    *       public void actions() {  \dots  }
    *
    * @endcode
    *
    *  and then invoke the constructor "<tt>new Bang (10.0)</tt>" to get a
    * `Bang` in 10 units of time. This is equivalent to "<tt>new
    * Bang().schedule(10.0)</tt>."
    *  @param delay        simulation time that must pass before the event
    *                      happens
    */
   public Event() {
      this (Simulator.getDefaultSimulator());
   }

   /**
    * Construct a new event instance associated with the given simulator.
    *  @param sim          Instance of class Simulator associated with the
    *                      new Event
    */
   public Event (Simulator sim) {
      if (sim == null)
         throw new NullPointerException();
      eventTime = -10.0;
      priority = 1.0;
      this.sim = sim;
   }

   /**
    * Schedules this event to happen in `delay` time units, i.e., at time
    * `sim.time() + delay`, by inserting it in the event list. When two or
    * more events are scheduled to happen at the same time and with the
    * same priority, they are placed in the event list (and executed) in
    * the same order as they have been scheduled. Note that the priority
    * of this event should be adjusted using  #setPriority *before* it is
    * scheduled.
    *  @param delay        simulation time that must pass before the event
    *                      happens
    */
   public void schedule (double delay) {
      if (delay < 0.0)
         throw new IllegalArgumentException ("Cannot schedule in the past.");
      if (eventTime > -1.0)
         throw new IllegalStateException ("Event already scheduled");
      eventTime = sim.time() + delay;
      sim.eventList.add (this);
   }

   /**
    * Schedules this event as the *first* event in the event list, to be
    * executed at the current time (as the next event).
    */
   public void scheduleNext() {
      if (eventTime > -1.0)
         throw new IllegalStateException ("Event already scheduled");
      eventTime = sim.time();
      priority  = 0.0;
      sim.eventList.addFirst (this);
   }

   /**
    * Schedules this event to happen just before, and at the same time, as
    * the event `other`. For example, if `Bing` and `Bang` are `Event`
    * subclasses, after the statements
    * @code
    *
    *          Bang bigOne = new Bang().schedule(12.0);
    *          new Bing().scheduleBefore(bigOne);
    *
    * @endcode
    *  the event list contains two new events scheduled to happen in 12
    * units of time: a `Bing` event, followed by a `Bang` called `bigOne`.
    *  @param other        event before which this event will be scheduled
    */
   public void scheduleBefore (Event other) {
      if (eventTime > -1.0)
         throw new IllegalStateException ("Event already scheduled");
      eventTime = other.eventTime;
      priority = other.priority;
      sim.eventList.addBefore (this, other);
   }

   /**
    * Schedules this event to happen just after, and at the same time, as
    * the event `other`.
    *  @param other        event after which this event will be scheduled
    */
   public void scheduleAfter (Event other) {
      if (eventTime > -1.0)
         throw new IllegalStateException ("Event already scheduled");
      eventTime = other.eventTime;
      priority = other.priority;
      sim.eventList.addAfter (this, other);
   }

   /**
    * Cancels this event and reschedules it to happen in `delay` time
    * units.
    *  @param delay        simulation time units that must elapse before
    *                      the event happens
    */
   public void reschedule (double delay) {
      if (delay < 0.0)
         throw new IllegalArgumentException ("Cannot schedule in the past.");
      if (eventTime < -1.0)
         throw new IllegalStateException ("Event not scheduled");
      sim.getEventList().remove (this);
      eventTime = sim.time() + delay;
      sim.getEventList().add (this);
   }

   /**
    * Cancels this event before it occurs. Returns `true` if cancellation
    * succeeds (this event was found in the event list), `false`
    * otherwise.
    *  @return `true` if the event could be cancelled
    */
   public boolean cancel() {
      boolean removed = false;
      if (eventTime >= sim.time()) removed = sim.getEventList().remove (this);
      eventTime = -10.0;
      return removed;
   }

   /**
    * Finds the first occurence of an event of class "type" in the event
    * list, and cancels it. Returns `true` if cancellation succeeds,
    * `false` otherwise.
    *  @param type         name of an event subclass
    *  @return `true` if an event of this class was found and cancelled
    */
   public final boolean cancel (String type) {
      Event ev = sim.getEventList().getFirstOfClass (type);
      return ev.cancel();
   }

   /**
    * Returns the simulator linked to this event.
    *  @return the simulator linked to the event
    */
   public final Simulator simulator() {
      return sim;
   }

   /**
    * Sets the simulator associated with this event to `sim`. This method
    * should not be called while this event is in an event list.
    *  @param sim          the Simulator
    */
   public final void setSimulator (Simulator sim) {
      if (sim == null)
          throw new NullPointerException();
      if (eventTime > -1.0)
         throw new UnsupportedOperationException (
            "Unable to set Simulator, current Event already scheduled");
      this.sim = sim;
   }

   /**
    * Returns the (planned) time of occurence of this event.
    *  @return the time of occurence of the event
    */
   public final double time() {
      return eventTime;
   }

   /**
    * Sets the (planned) time of occurence of this event to `time`. This
    * method should never be called after the event was scheduled,
    * otherwise the events would not execute in ascending time order
    * anymore.
    *  @param time         new time of occurence for the event
    */
   public final void setTime (double time) {
      if (eventTime > -1.0)
         throw new UnsupportedOperationException(
            "Unable to set time, current Event already scheduled");
      eventTime = time;
   }

   /**
    * Returns the priority of this event.
    *  @return the priority of the event
    */
   public final double priority() {
      return priority;
   }

   /**
    * Sets the priority of this event to `inPriority`. This method should
    * never be called after the event was scheduled, otherwise the events
    * would not execute in ascending priority order anymore.
    *  @param priority     new priority for the event
    */
   public final void setPriority (double priority) {
      if(eventTime > -1.0)
         throw new UnsupportedOperationException(
            "Unable to set priority, current Event already scheduled");
      this.priority = priority;
   }

   /**
    * Compares this object with the specified object `e` for order.
    * Returns @f$-1@f$ or @f$+1@f$ as this event occurs before or after
    * the specified event `e`, respectively. If the two events occur at
    * the same time, then returns @f$-1@f$, @f$0@f$, or @f$+1@f$ as this
    * event has a smaller, equal, or larger priority than event `e`.
    */
   public int compareTo (Event e) {
      if (eventTime < e.time())
         return -1;
      if (eventTime > e.time())
         return 1;
      // Si le moment de declenchement des "Event" est identique, on
      // examine leurs priorites.
      if (priority < e.priority())
         return -1;
      if (priority > e.priority())
         return 1;
      return 0;
   }

   /**
    * This is the method that is executed when this event occurs. Every
    * subclass of `Event` that is to be instantiated must provide an
    * implementation of this method.
    */
   public abstract void actions();

}