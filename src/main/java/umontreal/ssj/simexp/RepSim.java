/*
 * Class:        RepSim
 * Description:  simulation using independent runs or replications
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
 * @since        2007

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.simexp;

import umontreal.ssj.simevents.Simulator;

/**
 * Performs a simulation experiment on a finite horizon, using a certain
 * number of independent runs or replications. During each run&nbsp;@f$r@f$,
 * a complete simulation is executed, and the vector @f$\mathbf{X}_r@f$ is
 * generated. If simulation runs are independent, and the same system is
 * simulated during each run, after @f$n@f$ runs are performed, a sample
 * @f$(\mathbf{X}_0, …, \mathbf{X}_{n-1})@f$ is obtained.
 *
 * For such a simulation to be implemented, this class must be extended to
 * override the required methods:  #initReplicationProbes to initialize the
 * statistical probes collecting @f$\mathbf{X}_r’s@f$,  #initReplication(int)
 * to initialize the simulation model at the beginning of each replication,
 * and  #addReplicationObs(int) to add @f$\mathbf{X}_r@f$ to the statistical
 * probes.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class RepSim extends SimExp {
   private int minReps;
   private int maxReps;
   private int targetReps;
   private int doneReps;

   /**
    * Constructs a new replications-based simulator with a minimal number
    * of runs, `minReps`, and no maximal number of runs.
    *  @param minReps      the minimal number of replications.
    *  @exception IllegalArgumentException if the minimal number of
    * replications is smaller than 0.
    */
   public RepSim (int minReps) {
      this (minReps, Integer.MAX_VALUE);
   }

   /**
    * Constructs a new replications-based simulator with a minimal number
    * of runs `minReps`, and a maximal number of runs `maxReps`. This
    * maximum is used to avoid too long simulations when using sequential
    * sampling.
    *  @param minReps      the minimal number of replications.
    *  @param maxReps      the maximal number of replications.
    *  @exception IllegalArgumentException if the minimal or maximal
    * numbers of replications is smaller than 0, or if `minReps` is
    * greater than `maxReps`.
    */
   public RepSim (int minReps, int maxReps) {
      this (Simulator.getDefaultSimulator(), minReps, maxReps);
   }

   /**
    * Equivalent to the first constructor, with the given simulator `sim`.
    *  @param sim          the simulator attached to this object.
    *  @param minReps      the minimal number of replications.
    *  @exception IllegalArgumentException if the minimal number of
    * replications is smaller than 0.
    */
   public RepSim (Simulator sim, int minReps) {
      this (sim, minReps, Integer.MAX_VALUE);
   }

   /**
    * Equivalent to the second constructor, with the given simulator
    * `sim`.
    *  @param sim          the simulator attached to this object.
    *  @param minReps      the minimal number of replications.
    *  @param maxReps      the maximal number of replications.
    *  @exception IllegalArgumentException if the minimal or maximal
    * numbers of replications is smaller than 0, or if `minReps` is
    * greater than `maxReps`.
    */
   public RepSim (Simulator sim, int minReps, int maxReps) {
      super (sim);
      if (minReps <= 0)
         throw new IllegalArgumentException ("minReps <= 0");
      if (minReps > maxReps)
         throw new IllegalArgumentException ("minReps > maxReps");
      this.minReps = minReps;
      this.maxReps = maxReps;
      targetReps = minReps;
   }

   /**
    * Returns the minimal number of replications to be simulated before an
    * error check.
    *  @return the minimal number of replications.
    */
   public int getMinReplications() {
      return minReps;
   }

   /**
    * Sets the minimal number of replications required before an error
    * check to `minReps`. This also updates the maximal number of
    * replications if this maximum is smaller than the new minimum. This
    * will take effect only at the next call to  #simulate.
    *  @param minReps      the minimal number of replications.
    *  @exception IllegalArgumentException if the specified number of
    * replications is negative.
    */
   public void setMinReplications (int minReps) {
     if (minReps <= 0)
        throw new IllegalArgumentException ("minReps <= 0");
     this.minReps = minReps;
     if (maxReps < minReps)
        maxReps = minReps; 
   }

   /**
    * Returns the maximal number of replications to be simulated before an
    * error check. By default, this is set to
    * java.lang.Integer.MAX_VALUE, which is equivalent to infinity in
    * practice.
    *  @return the maximal number of replications.
    */
   public int getMaxReplications() {
      return maxReps;
   }

   /**
    * Sets the maximal number of replications required before an error
    * check to `maxReps`. This will take effect only at the next call to
    * #simulate.
    *  @param maxReps      the maximal number of replications.
    *  @exception IllegalArgumentException if the specified number of
    * replications is negative.
    */
   public void setMaxReplications (int maxReps) {
      if (maxReps < minReps)
         throw new IllegalArgumentException ("maxReps < minReps");
      this.maxReps = maxReps;
   }

   /**
    * Returns the actual target number of replications to be simulated
    * before an error check. By default, this is initialized to the
    * minimal number of replications, and is increased if new replications
    * are needed. However, it is not decreased by default, even upon a new
    * call to  #simulate.
    *  @return the target number of replications.
    */
   public int getTargetReplications() {
      return targetReps;
   }

   /**
    * Sets the target number of simulated replications before an error
    * check to `targetReps`. The value of `targetReps` must not be smaller
    * than the minimal number of replications returned by
    * #getMinReplications, or greater than the maximal number of
    * replications returned by  #getMaxReplications.
    *  @param targetReps   the target number of replications.
    *  @exception IllegalArgumentException if the new target number of
    * replications is smaller than the minimal number of replications.
    */
   public void setTargetReplications (int targetReps) {
      if (targetReps < minReps)
         throw new IllegalArgumentException
            ("Target number of replications too small");
      if (targetReps > maxReps)
         throw new IllegalArgumentException
            ("Target number of replications too large");
      this.targetReps = targetReps;
   }

   /**
    * Returns the total number of completed replications for the current
    * experiment.
    *  @return the number of completed replications.
    */
   public int getCompletedReplications() {
      return doneReps;
   }

   /**
    * Initializes any statistical collector used to collect values for
    * replications.
    */
   public abstract void initReplicationProbes();

   /**
    * Contains the necessary logic to perform the <tt>r</tt>th replication
    * of the simulation. By default, the method calls
    * umontreal.ssj.simevents.Sim.init to clear the event list, and uses
    * #initReplication(int) to initialize the model. It then calls
    * umontreal.ssj.simevents.Sim.start to start the simulation, calls
    * #replicationDone to increment the number of completed replications,
    * and  #addReplicationObs(int) to add observations to statistical
    * probes.
    *  @param r            the index of the replication.
    */
   public void performReplication (int r) {
      simulator().init();
      initReplication (r);
      simulator().start();
      replicationDone();
      addReplicationObs (r);
   }

   /**
    * Increments by one the number of completed replications. This is used
    * by  #performReplication(int).
    */
   protected void replicationDone() {
      ++doneReps;
   }

   /**
    * Initializes the simulation model for a new replication&nbsp;`r`.
    * This method should reset any counter and model state, and schedule
    * needed events. After the method returns, the model should be ready
    * for calling  umontreal.ssj.simevents.Sim.start. This method is
    * called just after the simulator is initialized.
    *  @param r            the index of the replication.
    */
   public abstract void initReplication (int r);

   /**
    * Adds statistical observations for the replication&nbsp;`r`. This
    * method is called just after the replication&nbsp;`r` is simulated.
    *  @param r            the index of the replication.
    */
   public abstract void addReplicationObs (int r);

   /**
    * Returns the approximate number of additional replications to meet an
    * experiment-specific stopping criterion. This is called after
    * #getTargetReplications replications are simulated. Since sequential
    * sampling is not used by default, the default implementation returns
    * 0, which stops the simulation after  #getTargetReplications
    * replications.
    *  @return the number of required additional replications.
    */
   public int getRequiredNewReplications() {
      return 0;
   }

   /**
    * Initializes this simulator for a new experiment. This method resets
    * the number of completed replications to 0, and calls
    * #initReplicationProbes to initialize statistical probes. This method
    * is called by  #simulate.
    */
   public void init() {
      if (simulating)
         throw new IllegalStateException ("Already simulating");
      doneReps = 0;
      if (targetReps < minReps)
         // minReps has been increased, so increase targetReps too
         targetReps = minReps;
      initReplicationProbes();
   }

   /**
    * Adjusts the target number of replications to simulate
    * `numNewReplications` additional replications. This method increases
    * the target number of replications by `numNewReplications`, and sets
    * the target number of replications to  #getMaxReplications if the new
    * target exceeds the maximal number of replications. This is called by
    * #simulate for sequential sampling.
    *  @param numNewReplications the number of additionnal replications
    *                            needed.
    */
   public void adjustTargetReplications (int numNewReplications) {
      if (numNewReplications < 0)
         throw new IllegalArgumentException ("numReplications < 0");
      if (numNewReplications == 0)
         return;
      targetReps = doneReps + numNewReplications;
      if (targetReps > maxReps)
         targetReps = maxReps;
   }

   /**
    * Simulates several independent simulation replications of a system.
    * When this method is called, the method  #init is called to
    * initialize the system, and  #getTargetReplications replications are
    * simulated by using  #performReplication(int). When the target number
    * of replications is simulated, the stopping condition is checked
    * using  #getRequiredNewReplications, and the target number of
    * replications is adjusted using  #adjustTargetReplications(int).
    * Additional replications are simulated until the method
    * #getRequiredNewReplications returns 0, or  #getMaxReplications
    * replications are simulated.
    */
   public void simulate() {
      init();
      simulating = true;
      try {
         while (doneReps < targetReps) {
            for (int i = 0; i < targetReps; i++)
               performReplication (i);
            adjustTargetReplications (getRequiredNewReplications());
         }
      }
      finally {
         simulating = false;
      }
   }


   public String toString() {
      StringBuffer sb = new StringBuffer (getClass().getName());
      sb.append ('[');
      sb.append ("minimal number of replications: ").append (minReps);
      if (maxReps < Integer.MAX_VALUE)
         sb.append (", maximal number of replications: ").append (maxReps);
      sb.append (", target number of replications: ").append (targetReps);
      if (simulating)
         sb.append (", simulation in progress");
      else
         sb.append (", simulation stopped");
      sb.append (", number of completed replications: ").append (doneReps);
      sb.append (']');
      return sb.toString();
   }
}