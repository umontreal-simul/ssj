/*
 * Class:        Accumulate
 * Description:  collects statistics on a variable that evolves in
                 simulation time
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
// This class doesn't belong to package stat because objects of this class
// always depend of Simulator

import java.util.Observable;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.stat.StatProbe;

/**
 * A subclass of  @ref umontreal.ssj.stat.StatProbe, for collecting
 * statistics on a variable that evolves in simulation time, with a
 * piecewise-constant trajectory. Each time the variable changes its value,
 * the method  #update(double) must be called to inform the probe of the new
 * value. The probe can be reinitialized by  #init.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Accumulate extends StatProbe implements Cloneable {

   private double initTime;    // Initialization time.
   private double lastTime;    // Last update time.
   private double lastValue;   // Value since last update.
   private Simulator sim;

   /**
    * Constructs a new `Accumulate` statistical probe using the default
    * simulator and initializes it by invoking `init()`.
    */
   public Accumulate() {
      super();
      sim = Simulator.getDefaultSimulator();
      init();
   }

   /**
    * Constructs a new `Accumulate` statistical probe linked to the given
    * simulator, and initializes it by invoking `init()`.
    *  @param inSim        the simulator of the current variable
    */
   public Accumulate (Simulator inSim) {
      super();
      if (inSim == null)
          throw new NullPointerException();
      sim = inSim;
      init();
   }

   /**
    * Constructs and initializes a new `Accumulate` statistical probe with
    * name `name` and initial time 0, using the default simulator.
    */
   public Accumulate (String name) {
      super();
      sim = Simulator.getDefaultSimulator();
      this.name = name;
      init();
   }

   /**
    * Constructs-initializes a new `Accumulate` statistical probe with
    * name `name` and initial time 0.
    *  @param name         descriptive name for the probe
    *  @param inSim        the simulator of the current variable
    */
   public Accumulate (Simulator inSim, String name) {
      super();
      if (inSim == null)
          throw new NullPointerException();
      sim = inSim;
      this.name = name;
      init();
   }

   /**
    * Initializes the statistical collector and puts the current value of
    * the corresponding variable to 0. **Note:** the initialization time,
    * the last update time and the simulation time are not reset to 0 by
    * this method. For this, `Sim.init()` must be used.
    */
   public void init() {
       maxValue = Double.MIN_VALUE;
       minValue = Double.MAX_VALUE;
       lastValue = 0.0;
       sumValue = 0.0;
       // May start the accumulator at t > 0; for ex., a warm-up period or
       // other reasons
       initTime = lastTime = sim.time();
   }

   /**
    * Same as  #init followed by  {@link #update(double) update(x)}.
    *  @param x            initial value of the probe
    */
   public void init (double x) {
       init();  update (x);
   }

   /**
    * Updates the accumulator using the last value passed to
    * #update(double).
    */
   public void update() {
      update (lastValue);
   }

   /**
    * Gives a new observation `x` to the statistical collector. If
    * broadcasting to observers is activated for this object, this method
    * will also transmit the new information to the registered observers
    * by invoking the methods  #notifyListeners(double).
    *  @param x            new observation given to the probe
    */
   public void update (double x) {
      if (collect) {
         double time = sim.time();
         if (x < minValue) minValue = x;
         if (x > maxValue) maxValue = x;
         sumValue += lastValue * (time - lastTime);
         lastValue = x;
         lastTime = time;
      }
      if (broadcast) {
         //setChanged();
         notifyListeners (x);
      }
   }
public double sum() {
      update (lastValue);
      return sumValue;
   }

/**
 * Returns the time-average since the last initialization to the last call to
 * `update`.
 */
public double average() {
      update (lastValue);
      double periode = lastTime - initTime;
      if (periode > 0.0)  return sumValue/periode;
      else  return 0.0;
   }

   public String shortReportHeader() {
      PrintfFormat pf = new PrintfFormat();
      pf.append (-9, "from time").append ("   ");
      pf.append (-9, "to time").append ("   ");
      pf.append (-8, "   min").append ("   ");
      pf.append (-8, "   max").append ("   ");
      pf.append (-10, " average");
      return pf.toString();
   }

   public String shortReport() {
      update();
      PrintfFormat pf = new PrintfFormat();
      pf.append (9, 2, 2, getInitTime()).append ("   ");
      pf.append (9, 2, 2, getLastTime()).append ("   ");
      pf.append (8, 3, 2, min()).append ("   ");
      pf.append (8, 3, 2, max()).append ("   ");
      pf.append (10, 3, 2, average());
      return pf.toString();
   }

   /**
    * Returns a string containing a report on this collector since its
    * last initialization.
    */
   public String report() {
      update (lastValue);
      PrintfFormat str = new PrintfFormat();
      str.append ("REPORT on Accumulate stat. collector ==> " + name);
      str.append (PrintfFormat.NEWLINE + "      from time   to time       min         max");
      str.append ("         average").append(PrintfFormat.NEWLINE);
      str.append (12, 2, 2, initTime);
      str.append (13, 2, 2, lastTime);
      str.append (11, 3, 2, minValue);
      str.append (12, 3, 2, (double)maxValue);
      str.append (14, 3, 2, (double)average()).append (PrintfFormat.NEWLINE);

      return str.toString();
    }

   /**
    * Returns the initialization time for this object. This is the
    * simulation time when  #init was called for the last time.
    *  @return the initialization time for this object
    */
   public double getInitTime() {
      return initTime;
   }

   /**
    * Returns the last update time for this object. This is the simulation
    * time of the last call to  #update or the initialization time if
    * #update was never called after  #init.
    *  @return the last update time of this object
    */
   public double getLastTime() {
      return lastTime;
   }

   /**
    * Returns the value passed to this probe by the last call to its
    * #update method (or the initial value if  #update was never called
    * after  #init ).
    *  @return the last update value for this object
    */
   public double getLastValue() {
      return lastValue;
   }

   /**
    * Returns the simulator associated with this statistical probe.
    *  @return the associated simulator.
    */
   public Simulator simulator() {
      return sim;
   }

   /**
    * Sets the simulator associated with this probe to `sim`. One should
    * call  #init after this method to reset the statistical probe.
    *  @param sim          the simulator of this probe
    */
   public void setSimulator(Simulator sim) {
       if (sim == null)
          throw new NullPointerException();
      this.sim = sim;
   }

   /**
    * Clone this object.
    */
   public Accumulate clone() {
      try {
         return (Accumulate)super.clone();
      } catch (CloneNotSupportedException e) {
         throw new IllegalStateException ("Accumulate can't clone");
      }
   }

}