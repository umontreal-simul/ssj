/*
 * Class:        Continuous
 * Description:  provides the basic structures and tools for 
                 continuous-time simulation
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
 * Represents a variable in a continuous-time simulation. This abstract class
 * provides the basic structures and tools for continuous-time simulation,
 * where certain variables evolve continuously in time, according to
 * differential equations. Such continuous variables can be mixed with events
 * and processes.
 *
 * Each type of continuous-time variable should be defined as a subclass of
 * `Continuous`. The instances of these subclasses are the actual
 * continuous-time variables. Each subclass must implement the method
 * #derivative which returns its derivative with respect to time. The
 * trajectory of this variable is determined by integrating this derivative.
 * The subclass may also reimplement the method  #afterEachStep, which is
 * executed immediately after each integration step.  By default (in the
 * class <tt>Continuous</tt>), this method does nothing. This method could,
 * for example, verify if the variable has reached a given threshold, or
 * update a graphical illustration of the variable trajectory.
 *
 * When creating a class representing a continuous variable, the  #toString
 * method can be overridden to display information about the continuous
 * variable. This information will be displayed when formating the event list
 * as a string.
 *
 * Each continuous variable has a linked simulator represented by an instance
 * of the  @ref Simulator class. If no simulator is provided explicitly when
 * constructing a variable, the default simulator returned by
 * `Simulator.getDefaultSimulator` is used.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class Continuous {

   // Private variables:

   boolean active; // This variable is currently being integrated.
   double  value;  // Current value of the variable.
   Event  ev;      // Event to be executed after each integ. step,

   //String name;
   double phi;
   double pi;
   double buffer;
   double sum;

   private Simulator sim;

   /**
    * Constructs a new continuous-time variable linked to the default
    * simulator, *without* initializing it.
    */
   public Continuous() {
      active = false;
      this.sim = Simulator.getDefaultSimulator();
   }

   /**
    * Constructs a new continuous-time variable linked to the given
    * simulator, *without* initializing it.
    *  @param sim          the simulator associated to this variable.
    *
    * @code
    *
    *    public Continuous (String name) \begin{hide} {
    *       this();
    *       this.name = name;
    *    } \end{hide}
    *
    * @endcode
    *  Constructs a new continuous-time variable (same as
    * <tt>Continuous()</tt>) with name `name` and linked to the default
    * simulator. This name can be used to identify the `Continuous`
    * variable in traces and reports.
    * @code
    *
    *    public Continuous (Simulator sim, String name) \begin{hide} {
    *       this(sim);
    *       this.name = name;
    *    } \end{hide}
    *
    * @endcode
    *  Constructs a new continuous-time variable (same as
    * <tt>Continuous(sim)</tt>) with name `name`. This name can be used to
    * identify the `Continuous` variable in traces and reports.
    */
   public Continuous (Simulator sim) {
       if (sim == null)
          throw new NullPointerException();
      active = false;
      this.sim = sim;
   }

   /**
    * Initializes or reinitializes the continuous-time variable to `val`.
    *  @param val          initial value of the variable
    */
   public void init (double val) {
      value = val;
   }

   /**
    * Returns the current value of this continuous-time variable.
    *  @return the current value of the variable
    */
   public double value() {
      return value;
   }

   /**
    * Returns the simulator linked to this continuous-time variable.
    *  @return the current simulator of the variable
    */
   public Simulator simulator() {
      return sim;
   }

   /**
    * Sets the simulator linked to this continuous-time variable. This
    * method should not be called while this variable is active.
    *  @param sim          the simulator of the current variable
    */
   public void setSimulator(Simulator sim) {
       if (sim == null)
          throw new NullPointerException();
      this.sim = sim;
   }

   /**
    * Starts the integration process that will change the state of this
    * variable at each integration step.
    */
   public void startInteg() {
      sim.continuousState().startInteg(this);
   }

   /**
    * Same as  #startInteg, after initializing the variable to `val`.
    *  @param val          initial value to start integration from
    */
   public void startInteg (double val) {
      init (val);   startInteg();
   }

   /**
    * Stops the integration process for this continuous variable. The
    * variable keeps the value it took at the last integration step before
    * calling `stopInteg`.
    */
   public void stopInteg() {
      sim.continuousState().stopInteg(this);
   }

   /**
    * This method should return the derivative of this variable with
    * respect to time, at time @f$t@f$. Every subclass of `Continuous`
    * that is to be instantiated must implement it. If the derivative does
    * not depend explicitly on time, @f$t@f$ becomes a dummy parameter.
    * Internally, the method is used with @f$t@f$ not necessarily equal to
    * the current simulation time.
    *  @param t            time at which the derivative must be computed
    */
   public abstract double derivative (double t);

   /**
    * This method is executed after each integration step for this
    * `Continuous` variable. Here, it does nothing, but every subclass of
    * `Continuous` may reimplement it.
    */
   public void afterEachStep() {
   }

   /**
    * Selects the Euler method as the integration method, with the
    * integration step size `h`, in time units, for the default simulator.
    * The non-static method  #selectEuler in  @ref ContinuousState can be
    * used to set the integration method for any given simulator. This
    * method appears here only to keep compatibility with older versions
    * of SSJ; using a non-static  @ref Simulator instance rather than the
    * default simulator is recommended.
    *  @param h            integration step, in simulation time units
    */
   public static void selectEuler(double h) {
      Simulator.getDefaultSimulator().continuousState().selectEuler(h);
   }

   /**
    * Selects a Runge-Kutta method of order&nbsp;4 as the integration
    * method to be used, with step size `h`. The non-static method
    * #selectRungeKutta4 in  @ref ContinuousState can be used to set the
    * integration method for any given simulator. This method appears here
    * only to keep compatibility with older versions of SSJ; using a
    * non-static  @ref Simulator instance rather than the default
    * simulator is recommended.
    */
   public static void selectRungeKutta4(double h) {
      Simulator.getDefaultSimulator().continuousState().selectRungeKutta4(h);
   }

   /**
    * Selects a Runge-Kutta method of order&nbsp;2 as the integration
    * method to be used, with step size `h`. The non-static method
    * #selectRungeKutta2 in  @ref ContinuousState can be used to set the
    * integration method for any given simulator. This method appears here
    * only to keep compatibility with older versions of SSJ; using a
    * non-static  @ref Simulator instance rather than the default
    * simulator is recommended.
    */
   public static void selectRungeKutta2(double h) {
      Simulator.getDefaultSimulator().continuousState().selectRungeKutta2(h);
   }

}