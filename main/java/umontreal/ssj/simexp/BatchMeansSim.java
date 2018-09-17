/*
 * Class:        BatchMeansSim
 * Description:  
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

import java.lang.ref.SoftReference;

import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

import umontreal.ssj.simevents.Event;
import umontreal.ssj.simevents.Simulator;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.util.Misc;

/**
 * Performs a simulation experiment on an infinite horizon, for estimating
 * steady-state performance measures, using batch means. Batches are
 * delimited using a user-specified condition such as a fixed duration in
 * simulation time units (the default), the number of occurrences of an event
 * such as the arrival of a customer, a regenerative cycle, etc. After the
 * condition for batch termination is defined, the *batch size* can be set.
 * This size can be, depending on how batches are delimited, a time duration,
 * a number of events, or 1 for regenerative cycles. The *batch length* is
 * defined to be the duration of a batch, in simulation time units,
 * independently of how batches are defined. By default, the batch size and
 * the batch length are equivalent (and constant), but they may differ if the
 * condition for batch termination is changed.
 *
 * A warmup period is usually simulated to reduce the bias induced by the
 * initial state of the system. During the warmup, the system runs without
 * any observation being collected. By default, this period has a duration
 * fixed in simulation time units, but this can also be changed.
 *
 * After the warmup is over, events are counted as follows. The simulation
 * model defines some counters being updated when events occur and reset at
 * the beginning of each batch. At the end of a batch, the value of these
 * counters are used to generate a random vector @f$\mathbf{V}_j@f$ before
 * the counters are reset. Alternatively, a simulation may compute and update
 * the @f$\mathbf{V}_j@f$’s directly, without using intermediate counters. A
 * set of data structures is needed to collect and store these
 * @f$\mathbf{V}_j@f$’s, the simplest option being a set of
 * @ref umontreal.ssj.stat.TallyStore instances. This generates values for
 * @f$m@f$ *real batches*. The sample @f$(\mathbf{X}_0,…,\mathbf{X}_{n-1})@f$
 * is then obtained from @f$(\mathbf{V}_0,…,\mathbf{V}_{m-1})@f$, so a second
 * set of data structures is needed to collect the @f$\mathbf{X}_r@f$’s, the
 * simplest being a set of  @ref umontreal.ssj.stat.Tally instances. The most
 * straightforward way to estimate covariances on components of
 * @f$\mathbf{X}_r@f$ is by considering the vectors @f$\mathbf{X}_r@f$
 * i.i.d., which is not true in general. However, by choosing a sufficiently
 * large simulation time and batch lengths, the correlation between batches
 * can be reduced, and the @f$\mathbf{X}_r@f$’s are approximately i.i.d. and
 * normally distributed. If batches correspond to regenerative cycles, the
 * @f$\mathbf{X}_r@f$ are then truly i.i.d., but only approximately normally
 * distributed. Confidence intervals on functions of @f$\boldsymbol{\mu}@f$
 * can be approximated using the central limit theorem, as with independent
 * replications, or using the delta theorem for functions of multiple means
 * or when batches have different lengths.
 *
 * The sample size corresponding to the number of simulated batches is always
 * fixed when sequential sampling is not used; we have @f$n=m@f$ and
 * @f$\mathbf{X}_r=\mathbf{V}_r@f$ for @f$r=0,…,n-1@f$. However, when
 * sequential sampling is used, the sample size can be random or fixed. If
 * the sample size is random, we still have @f$\mathbf{X}_r=\mathbf{V}_r@f$.
 * When using this mode, the batch size must be chosen carefully to reduce
 * correlation between batches.
 *
 * On the other hand, if the sample size @f$n@f$ is required to be fixed
 * while sequential sampling is used, *batch aggregation* must be activated
 * to have *effective batches* with random lengths. When aggregation is
 * enabled, real batches are simulated as usual, but the obtained values
 * @f$\mathbf{V}_j@f$ are regrouped (or aggregated) to form effective
 * batches. The number of real batches must be @f$m=h*n@f$ to get a sample of
 * size @f$n@f$, and @f$h@f$ can be any integer greater than or equal to 1.
 * In this case,
 * @f[
 *   \mathbf{X}_r=\sum_{k=0}^{h-1} \mathbf{V}_{hr+k}.
 * @f]
 * The size of the effective batches increases while the sample size remains
 * fixed because the group size @f$h@f$ increases with simulation length.
 * When aggregation is used, the effective batches become longer while
 * simulation time increases, and the correlation between effective batches
 * should decrease with simulation time. However, the increment to the target
 * number of batches must always be a multiple of the sample size, which can
 * lead to useless simulation if the sample size is large.
 *
 * This class must be extended to implement a batch means simulator, and the
 * appropriate methods must be defined or overridden. This class uses a
 * simulation event to stop the simulation at the end of the warmup period
 * and batches. One can use this event for fixed-duration warmup and batches,
 * or schedule their own events which call `Sim.stop` to end warmup or
 * batches. To change how the warmup period is terminated, one must override
 * the method  #warmup. For the batch termination condition to be redefined,
 * #simulateBatch must be overridden.
 *
 * One must implement  #initSimulation to initialize the simulated model
 * before the warmup,  #initBatchStat to reset the model-specific counters
 * used to compute the @f$\mathbf{V}_j@f$’s,  #initRealBatchProbes and
 * #addRealBatchObs to initialize statistical probes and add observations for
 * real batches,  #initEffectiveBatchProbes and
 * #addEffectiveBatchObs(int,int,double) to initialize statistical probes and
 * add observations for effective batches.
 *
 * The moment the latter methods are called depends on the status of batch
 * aggregation: when aggregation is turned ON, before any error check or the
 * end of the simulation,  #initEffectiveBatchProbes is called once before
 * #addEffectiveBatchObs(int,int,double) is called @f$n@f$ successive times,
 * with different parameters. When aggregation is turned OFF, then method
 * #initEffectiveBatchProbes is called after the warmup is over, and
 * #addEffectiveBatchObs(int,int,double) is called each time a batch ends,
 * after  #addRealBatchObs is called.
 *
 * If sequential sampling is used,  #getRequiredNewBatches must be overridden
 * to implement error checking. In some particular situations, the user may
 * also need to override  #allocateCapacity(int) and
 * #regroupRealBatches(int).
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class BatchMeansSim extends SimExp {
   private DoubleArrayList batchTimes = new DoubleArrayList();
   private boolean batchLengths;
   private boolean aggregation;
   private int minBatches;
   private int maxBatches;
   private double warmupTime;
   private EndSimEvent endSimEvent;

   private boolean warmupDone;
   private int targetBatches;
   private int doneBatches;
   private int droppedBatches;
   private double batchFraction;
   private double batchSizeMultiplier;
   private boolean aggregateUpdated;
   private double batchSize;

   private int nAgr;

   /**
    * Constructs a new batch means simulator using at least `minBatches`
    * batches with size `batchSize`, with a warmup period of duration
    * `warmupTime`. By default, batch aggregation and batch lengths
    * keeping are turned off, and the maximal number of batches is
    * infinite.
    *  @param minBatches   the minimal number of batches to simulate.
    *  @param batchSize    the size of the batches.
    *  @param warmupTime   the duration of the warmup period.
    *  @exception IllegalArgumentException if the minimal number of
    * batches is smaller than or equal to 0, or if the warmup time is
    * smaller than 0.
    */
   public BatchMeansSim (int minBatches, double batchSize,
                         double warmupTime) {
      this (minBatches, Integer.MAX_VALUE, batchSize, warmupTime);
   }

   /**
    * Constructs a batch means simulator with a maximum of `maxBatches`
    * batches to avoid excessive memory usage and too long simulations
    * when using sequential sampling. See
    * #BatchMeansSim(int,double,double) for more information about the
    * other parameters.
    *  @param minBatches   the minimal number of batches to simulate.
    *  @param maxBatches   the maximal number of batches to simulate.
    *  @param batchSize    the size of the batches.
    *  @param warmupTime   the duration of the warmup period.
    *  @exception IllegalArgumentException if the minimal number of
    * batches is smaller than or equal to 0, or if the warmup time is
    * smaller than 0.
    */
   public BatchMeansSim (int minBatches, int maxBatches, double batchSize,
                         double warmupTime) {
      this (Simulator.getDefaultSimulator(), minBatches, maxBatches,
      batchSize, warmupTime);
   }

   /**
    * Equivalent to the first constructor, with a user-defined simulator
    * `sim`.
    *  @param sim          the simulator attached to this object.
    *  @param minBatches   the minimal number of batches to simulate.
    *  @param batchSize    the size of the batches.
    *  @param warmupTime   the duration of the warmup period.
    *  @exception IllegalArgumentException if the minimal number of
    * batches is smaller than or equal to 0, or if the warmup time is
    * smaller than 0.
    */
   public BatchMeansSim (Simulator sim, int minBatches, double batchSize,
                         double warmupTime) {
      this (sim, minBatches, Integer.MAX_VALUE, batchSize, warmupTime);
   }

   /**
    * Equivalent to the second constructor, with a user-defined simulator
    * `sim`.
    *  @param sim          the simulator attached to this object.
    *  @param minBatches   the minimal number of batches to simulate.
    *  @param maxBatches   the maximal number of batches to simulate.
    *  @param batchSize    the size of the batches.
    *  @param warmupTime   the duration of the warmup period.
    *  @exception IllegalArgumentException if the minimal number of
    * batches is smaller than or equal to 0, or if the warmup time is
    * smaller than 0.
    */
   public BatchMeansSim (Simulator sim, int minBatches, int maxBatches,
                         double batchSize, double warmupTime) {
      super (sim);
      if (minBatches <= 0)
         throw new IllegalArgumentException
            ("minBatches <= 0");
      if (maxBatches < minBatches)
         throw new IllegalArgumentException ("maxBatches < minBatches");
      if (warmupTime < 0)
         throw new IllegalArgumentException
            ("Warmup time must not be negative");
      if (batchSize <= 0)
         throw new IllegalArgumentException
            ("Batch size must not be 0 or negative");
      this.minBatches = minBatches;
      this.maxBatches = maxBatches;
      this.warmupTime = warmupTime;
      this.batchSize = batchSize;
      targetBatches = minBatches;
   }

   /**
    * Returns `true` if the aggregation of batches is turned ON. If
    * #getRequiredNewBatches always returns 0 (the default), the
    * aggregation has no effect since the number of batches is not random.
    * By default, batch aggregation is turned OFF.
    *  @return the batch aggregation indicator.
    */
   public boolean getBatchAggregation() {
      return aggregation;
   }

   /**
    * Sets the batch aggregation indicator to `a`. This should not be
    * called during an experiment.
    *  @param a            the new batch aggregation indicator.
    *  @exception IllegalStateException if the warmup period is over, and
    * simulation is not terminated.
    */
   public void setBatchAggregation (boolean a) {
      if (warmupDone && simulating)
         throw new IllegalStateException
            ("Cannot change the aggregation status during simulation");
      aggregation = a;
   }

   /**
    * Indicates that the length, in simulation time units, of each real
    * batch has to be kept. By default, this is set to `false`. When batch
    * aggregation is turned ON, the batch lengths are always kept.
    *  @return the batch lengths keeping indicator.
    */
   public boolean getBatchLengthsKeeping() {
      return batchLengths;
   }

   /**
    * Sets the batch lengths keeping indicator to `b`. This has no impact
    * if batch aggregation is turned ON.
    *  @param b            the new value of the indicator.
    *  @exception IllegalStateException if the warmup period is over and
    * simulation is not terminated.
    */
   public void setBatchLengthsKeeping (boolean b) {
      if (warmupDone && simulating)
         throw new IllegalStateException
            ("Cannot change the aggregation status during simulation");
      batchLengths = b;
   }

   /**
    * Returns the minimal number of batches required for estimating the
    * steady-state performance measures of interest. If aggregation is
    * turned ON, this is also the final number of effective batches, i.e.,
    * the sample size.
    *  @return the minimal number of batches.
    */
   public int getMinBatches() {
      return minBatches;
   }

   /**
    * Sets the minimal number of batches to `minBatches`.
    *  @param minBatches   the new minimal number of batches.
    *  @exception IllegalArgumentException if the specified minimal number
    * of batches is smaller than 0.
    */
   public void setMinBatches (int minBatches) {
      if (minBatches <= 0)
         throw new IllegalArgumentException ("minBatches <= 0");
      this.minBatches = minBatches;
      if (maxBatches < minBatches)
         maxBatches = minBatches;
   }

   /**
    * Returns @f$M@f$, the maximal number of batches to be used for
    * estimating the steady-state performance measures of interest. This
    * is used to prevent the number of batches from growing indefinitely
    * when using sequential sampling. By default, this is set to
    * java.lang.Integer.MAX_VALUE, which is equivalent to infinity in
    * practice.
    *  @return the maximal number of batches.
    */
   public int getMaxBatches() {
      return maxBatches;
   }

   /**
    * Sets the maximal number of batches to `maxBatches`.
    *  @param maxBatches   the new maximal number of batches.
    *  @exception IllegalArgumentException if the maximal number of
    * batches is smaller than the actual minimal number of batches.
    */
   public void setMaxBatches (int maxBatches) {
      if (maxBatches < minBatches)
         throw new IllegalArgumentException ("maxBatches < minBatches");
      this.maxBatches = maxBatches;
   }

   /**
    * Returns the current batch size as defined for this simulator. By
    * default, this is a duration in simulation time units. Depending on
    * the batch termination condition, which can be changed by overriding
    * #simulateBatch, it can be the number of occurrences of an event, or
    * 1 for regenerative cycles.
    *  @return the current batch size.
    */
   public double getBatchSize() {
      return batchSize;
   }

   /**
    * Sets the batch size to `batchSize`.
    *  @param batchSize    the new batch size.
    *  @exception IllegalArgumentException if the given batch size is
    * negative or 0.
    */
   public void setBatchSize (double batchSize) {
      if (batchSize <= 0)
         throw new IllegalArgumentException
            ("batchSize <= 0");
      this.batchSize = batchSize;
   }

   /**
    * Returns the duration of the warmup period for the simulation. By
    * default, this duration is expressed in simulation time units, but
    * this can be changed by overriding  #warmup.
    *  @return the warmup time for the simulation.
    */
   public double getWarmupTime() {
      return warmupTime;
   }

   /**
    * Sets the warmup time to `warmupTime`. If this method is called while
    * a simulation is in progress, the new time will affect the next
    * simulation only.
    *  @param warmupTime   the new warmup time.
    *  @exception IllegalArgumentException if the warmup time is smaller
    * than 0.
    */
   public void setWarmupTime (double warmupTime) {
      if (warmupTime < 0)
         throw new IllegalArgumentException ("warmupTime < 0");
      this.warmupTime = warmupTime;
   }

   /**
    * Returns the remaining fraction of batch to be simulated. This method
    * is called when scheduling the end of the next batch during a
    * simulation. Sometimes, it can be necessary to increase the batch
    * size to avoid excessive memory usage. In this case, stored real
    * batches (the @f$\mathbf{V}_j@f$’s) have to be regrouped to use less
    * memory, and the last stored @f$\mathbf{V}_j@f$ may represent an
    * incomplete batch with respect to the new batch size. This method
    * returns the fraction of batch, with respect to the new batch size,
    * remaining to be simulated before a new batch starts. If regrouping
    * is not used, this always returns 1. The returned value is always
    * greater than 0 and smaller than or equal to 1. For example, with a
    * fixed time batch size `s`, the next end of batch will be scheduled
    * in `s*getBatchFraction()` simulation time units. This method returns
    * values different from 1 only if one overrides
    * #allocateCapacity(int), and  #regroupRealBatches.
    *  @return the remaining fraction of batch to simulate.
    *
    *  @see #regroupRealBatches(int)
    */
   public double getBatchFraction() {
      return batchFraction;
   }

   /**
    * Returns the batch size multiplier after the simulation of a new
    * batch. This is called when scheduling the end of a new batch, to
    * multiply the batch size after  #regroupRealBatches(int) has
    * regrouped real batches. This can return any value greater than 0, or
    * 1 if the size is unchanged (the most common case) or if aggregation
    * is not used. This method returns values different from 1 only if one
    * overrides  #allocateCapacity(int), and  #regroupRealBatches.
    *  @return the batch size multiplier.
    *
    *  @see #regroupRealBatches(int)
    */
   public double getBatchSizeMultiplier() {
      return batchSizeMultiplier;
   }

   /**
    * Returns the target number of simulated real batches at the next time
    * the stopping condition is checked. By default, this number is set to
    * the minimal number of batches and is increased if the stopping
    * condition check requires new simulated batches. This target number
    * of batches is not decreased automatically, even upon a new call to
    * #simulate.
    *  @return the target number of batches.
    */
   public int getTargetBatches() {
      return targetBatches;
   }

   /**
    * Sets the target number of simulated batches before an error check or
    * the end of the simulation to `targetBatches`.
    *  @param targetBatches the target number of batches.
    *  @exception IllegalArgumentException if the new target number of
    * batches is smaller than the minimal number of batches, or greater
    * than the maximal number of batches.
    */
   public void setTargetBatches (int targetBatches) {
      if (targetBatches < minBatches)
         throw new IllegalArgumentException
            ("Target number of batches too small");
      if (targetBatches > maxBatches)
         throw new IllegalArgumentException
            ("Target number of batches too large");
      this.targetBatches = targetBatches;
   }

   /**
    * Returns the number of completed real batches since the beginning of
    * the run.
    *  @return the number of completed real batches.
    */
   public int getCompletedRealBatches() {
      return doneBatches;
   }

   /**
    * Returns the number of real batches dropped. When using sequential
    * sampling, the target number of batches can become very high,
    * resulting in not enough memory available to store the real batches.
    * One simple heuristic to address this issue is to drop the first real
    * batches, and increase the batch length. This method gives the number
    * of real batches which have been dropped to save memory.
    *  @return the number of dropped real batches.
    */
   public int getDroppedRealBatches() {
      return droppedBatches;
   }

   /**
    * Drops the `n` first real batches to save memory.
    *  @param n            the number of real batches to drop.
    *  @exception IllegalArgumentException if `n` is negative or greater
    * than the number of stored real batches.
    *  @see #getDroppedRealBatches()
    */
   public void dropFirstRealBatches (int n) {
      if (n < 0 || n > doneBatches - droppedBatches)
         throw new IllegalArgumentException
         ("Cannot drop less than 0 or more than "
               + (doneBatches - droppedBatches) +
               " real batches");
      if (batchLengths || aggregation)
         batchTimes.removeFromTo (0, n - 1);
      droppedBatches += n;
   }

   /**
    * Returns the real batch corresponding to simulation time `time` when
    * batch lengths are kept. If batch lengths are not kept, or if the
    * given time corresponds to the warmup period, this method returns a
    * negative value.
    *  @param time         the simulation time.
    *  @return the batch corresponding to the time, or a negative value.
    */
   public int getBatch (double time) {
      if (!warmupDone)
         return -1;
      if (!batchLengths && !aggregation)
         return -1;
      if (batchTimes.size() == 0)
         return -1;
      return Misc.getTimeInterval (batchTimes.elements(), 0, doneBatches, time);
   }

   /**
    * Determines if the warmup period for the simulation is over.
    *  @return `true` if the warmup period is finished, `false` otherwise.
    */
   public boolean isWarmupDone() {
      return warmupDone;
   }

   /**
    * Returns @f$h@f$, the number of real batches contained into an
    * effective batch. If aggregation is turned OFF, this always returns 1
    * as soon as at least one batch is simulated. Otherwise, this returns
    * a number greater than or equal to 1 as soon as
    * #addEffectiveBatchObs(int,int,double) is called.
    *  @return the number of real batches in one effective batch.
    *
    *  @exception IllegalStateException if @f$h@f$ is not available.
    */
   public int getNumAggregates() {
      if (nAgr == 0)
         throw new IllegalStateException
            ("Number of real batches not available");
      return nAgr;
   }

   /**
    * Returns the length, in simulation time units, of the real batch
    * `batch`. If batch lengths are not kept, this method can return the
    * length of the last batch only.
    *  @param batch        the batch index.
    *  @return the length of the batch.
    *
    *  @exception IndexOutOfBoundsException if the batch index is out of
    * bounds.
    *  @exception IllegalArgumentException if the batch length is not
    * available for `batch`.
    */
   public double getRealBatchLength (int batch) {
      if (batchLengths || aggregation)
         // The first get will fail if batch + 1 is out of bounds.
         // If the first get succeeds, the second will, so we can
         // use getQuick to eliminate a check.
        return batchTimes.get (batch + 1 - droppedBatches) - batchTimes.getQuick (batch - droppedBatches);
      else if (batch != doneBatches - 1)
         throw new IllegalArgumentException ("Unavailable batch length");
      else
         return batchTimes.getQuick (1) - batchTimes.getQuick (0);
   }

   /**
    * Returns the starting simulation time of batch `batch`.
    *  @param batch        the queried batch index.
    *  @return the starting time.
    */
   public double getRealBatchStartingTime (int batch) {
      if (batchLengths || aggregation)
         return batchTimes.get (batch);
      else if (batch != doneBatches - 1)
         throw new IllegalArgumentException ("Unavailable batch time");
      else
         return batchTimes.getQuick (0);
   }

   /**
    * Returns the ending simulation time of batch `batch`.
    *  @param batch        the queried batch index.
    *  @return the ending time.
    */
   public double getRealBatchEndingTime (int batch) {
      if (batchLengths || aggregation)
         return batchTimes.get (batch + 1);
      else if (batch != doneBatches - 1)
         throw new IllegalArgumentException ("Unavailable batch time");
      else
         return batchTimes.getQuick (1);
   }

   /**
    * Allocates the necessary memory for storing `capacity` real batches.
    * When using sequential sampling, if the variance of estimators is
    * high, many additional batches may be needed to reach the target
    * precision. To avoid memory problems after a long simulation time,
    * this method can preallocate the necessary memory. The method must
    * ensure that the data structures used to store the
    * @f$\mathbf{V}_j@f$’s can contain `capacity` real batches. This is
    * done by recreating arrays or resizing data structures. By default,
    * this method throws an  UnsupportedOperationException.
    * #regroupRealBatches must be implemented if this method does not
    * throw this exception.
    *  @param capacity     the number of real batches to store.
    *  @exception UnsupportedOperationException to indicate that capacity
    * allocation is not supported
    */
   public void allocateCapacity (int capacity) {
      throw new UnsupportedOperationException();
   }

   /**
    * Regroups real batches `x` by `x`. When memory is low, the simulator
    * can try to regroup real batches and increase the batch size
    * consequently. This is partly done by this method. The user must
    * override it and modify the internal data structures storing the
    * @f$\mathbf{V}_j@f$’s in order to regroup elements. The number of
    * real batches @f$m@f$ becomes @f$m’=\lfloor m/x\rfloor@f$. After
    * this method returns, each new @f$\mathbf{V}_j@f$ should contain
    * @f$\mathbf{V}_j=\sum_{l=0}^{x-1}\mathbf{V}_{jx+l}@f$ for @f$j=0,…,
    * m’ - 1@f$, @f$\mathbf{V}_{m’}=\sum_{l=0}^{(m \bmod x) - 1}
    * \mathbf{V}_{m’x + l}@f$, and @f$\mathbf{V}_j=0@f$ for @f$j=m’+1, …,
    * m - 1@f$. Some static methods called  #regroupElements(double[],int)
    * are provided by this class to help the user with this. By default,
    * this method throws an  UnsupportedOperationException, disabling this
    * functionality which is not always needed.
    *  @param x            the regrouping factor.
    *  @exception UnsupportedOperationException if regrouping is not
    * supported.
    */
   public void regroupRealBatches (int x) {
      throw new UnsupportedOperationException();
   }


   private final void allocateCapacity() {
      // Preallocate the necessary memory to store
      // all the observations when using fixed
      // number of batches.
      // A maximal number of batches is used
      // to prevent memory overflow problems.
      // Unfortunately, in Java, we can only
      // use heuristics to prevent the memory
      // problem. In addition to a maximum
      // number of batches, we use a soft reference
      // that will be cleared before the OutOfMemoryError
      // is thrown.
      // The size of the referent can be tuned; it is only
      // an heuristic.
      SoftReference<?> softRef = new SoftReference<double[]> (new double[50000]);
      batchFraction = 1.0;
      batchSizeMultiplier = 1.0;
      int currentCapacity = doneBatches;
      while (currentCapacity < targetBatches && softRef.get() != null) {
         // Try to double the capacity
         int newCapacity = Math.min (2*currentCapacity + 1, targetBatches);
         // Since we want all collectors to have the same capacity
         // we must not stop the for loop if the softRef
         // is cleared.
         if (aggregation || batchLengths)
            batchTimes.ensureCapacity (newCapacity + 1);
         try {
            allocateCapacity (newCapacity);
         }
         catch (UnsupportedOperationException use) {
            return;
         }
         currentCapacity = newCapacity;
      }
      if (currentCapacity < targetBatches && softRef.get() == null) {
         if (!aggregation)
            throw new IllegalStateException
               ("Cannot increase batch size");
         // The soft reference was cleared.
         // The memory taken by the simulation
         // could crash the program.
         // To try avoiding that, we will stop taking memory
         // by increasing the batch size.
         double ftargetBatches = targetBatches;
         double fdoneBatches = doneBatches;
         while (targetBatches > currentCapacity) {
            batchSizeMultiplier *= 2;
            // From now, we can have a fraction of batches done
            fdoneBatches /= 2;
            ftargetBatches /= 2;
            targetBatches = (int)ftargetBatches;
            doneBatches = (int)fdoneBatches;
            if (aggregation || batchLengths) {
               for (int i = 0; i < batchTimes.size() / 2; i++)
                  batchTimes.setQuick (i, batchTimes.getQuick (2*i));
               batchTimes.setSize (batchTimes.size() / 2);
            }
            regroupRealBatches (2);
            // Compute the remaining fraction of batch to simulate
            batchFraction = 1.0 - (fdoneBatches - doneBatches);
         }
         // Round the new total number of batches to a multiple
         // of minBatches.
         if (targetBatches % minBatches != 0) {
            targetBatches /= minBatches;
            targetBatches *= minBatches;
            targetBatches += minBatches;
         }
      }
   }

/**
 * Initializes the simulator for a new run. This is called by the  #init
 * method after  umontreal.ssj.simevents.Sim.init is called.
 */
public abstract void initSimulation();

   /**
    * Resets the counters used for computing observations during the
    * simulation at the beginning of a new batch.
    */
   public abstract void initBatchStat();

   /**
    * Initializes any statistical collector for real batches. This is
    * called at the end of the warmup period.
    */
   public abstract void initRealBatchProbes();

   /**
    * Initializes any statistical collector for effective batches. This is
    * called at every stopping condition check when aggregation is ON, or
    * at the end of the warmup period when it is OFF.
    */
   public abstract void initEffectiveBatchProbes();

   /**
    * Collects values of a @f$\mathbf{V}_j@f$ vector concerning the last
    * simulated real batch. This method is called at the end of each real
    * batch.
    */
   public abstract void addRealBatchObs();

   /**
    * Adds an observation to each statistical collector corresponding to
    * an effective batch. The effective batch for which this method is
    * called has length `l`, and regroups real batches `s`, …, `s + h -
    * 1`. This method is called after each error check if aggregation is
    * turned ON, or after each real batch if it is turned OFF.
    */
   public abstract void addEffectiveBatchObs (int s, int h, double l);

   /**
    * Computes the approximate number of required real batches to be
    * simulated before the simulation can be stopped. The default
    * implementation always returns 0, which stops the simulation after
    * #getTargetBatches real batches are obtained; sequential sampling is
    * not used by default.
    *
    * Note: if the method uses
    * umontreal.ssj.simexp.SimExp.getRequiredNewObservations(StatProbe[],double,double)
    * with a statistical probe containing one observation per effective
    * batch, this gives the number of additional effective batches to
    * simulate. This value should be multiplied with  #getNumAggregates to
    * get the number of additional real batches.
    *  @return the approximate required number of additionnal real
    * batches.
    */
   public int getRequiredNewBatches() {
      return 0;
   }

   /**
    * Initializes the simulator for a new experiment. This method, called
    * by  #simulate, resets the counter for the number of batches, calls
    * `simulator().init`, followed by  #initSimulation.
    */
   public void init() {
      if (simulating)
         throw new IllegalStateException
         ("Already simulating");
      nAgr = 0;
      warmupDone = false;
      doneBatches = 0;
      droppedBatches = 0;
      if (targetBatches < minBatches)
         targetBatches = minBatches;
      if (aggregation && targetBatches % minBatches != 0) {
         targetBatches /= minBatches;
         targetBatches += minBatches;
         if (targetBatches + minBatches <= maxBatches)
            targetBatches += minBatches;
      }
      aggregateUpdated = false;
      simulator().init();
      endSimEvent = new EndSimEvent (simulator());
      initSimulation();
   }

   /**
    * Returns the event used to stop the simulation at the end of the
    * warmup or batches.
    */
   public Event getEndSimEvent() {
      return endSimEvent;
   }

   /**
    * Performs a warmup by calling  #warmup(double). By default, this
    * method calls  #warmup(double) with the value returned by
    * #getWarmupTime, but one can override this method to simulate the
    * warmup differently. If the duration of the warmup period is not
    * fixed, one can call  #warmup(double) with
    * <tt>Double.POSITIVE_INFINITY</tt>; this prevents the method from
    * scheduling the ending event, and let the simulator call
    * `simulator().stop` at appropriate time.
    */
   public void warmup() {
      warmup (warmupTime);
   }

   /**
    * Performs a warmup of fixed duration `warmupTime`. This method
    * simulates for `warmupTime` simulation time units, and initializes
    * statistical probes for real batches through  #initRealBatchProbes.
    * If the duration of the warmup period is not fixed, one can call
    * #warmup(double) with <tt>Double.POSITIVE_INFINITY</tt>; this
    * prevents the method from scheduling the ending event, and let the
    * simulator call `simulator().stop` at appropriate time.
    */
   public void warmup (double warmupTime) {
      if (warmupDone)
         throw new IllegalStateException ("Warmup already done");
      if (!Double.isInfinite (warmupTime) && !Double.isNaN (warmupTime))
         endSimEvent.schedule (warmupTime);
      simulator().start();
      endSimEvent.cancel();
      initRealBatchProbes();
      if (!aggregation)
         initEffectiveBatchProbes();
      if (aggregation || batchLengths) {
         batchTimes.setSize (0);
         batchTimes.add (simulator().time());
      }
      else {
         batchTimes.setSize (2);
         batchTimes.trimToSize();
         batchTimes.setQuick (0, simulator().time());
         batchTimes.setQuick (1, simulator().time());
      }
      warmupDone = true;
   }

   /**
    * Simulate a new batch with default length. By default, this method
    * multiplies the current batch size with  #getBatchSizeMultiplier, and
    * schedules the next end-batch event to happen in  #getBatchSize `*`
    * #getBatchFraction simulation time units, by using
    * #simulateBatch(double). After the batch is simulated, the batch-size
    * multiplier and fraction are reset to 1. If the batch lengths are not
    * fixed, one can override this method to call  #simulateBatch(double)
    * with <tt>Double.POSITIVE_INFINITY</tt>; this prevents the method
    * from scheduling the ending event, and let the simulator call
    * `simulator().stop` at appropriate time.
    */
   public void simulateBatch() {
      double f = getBatchSizeMultiplier();
      if (f > 1)
         batchSize *= f;
      simulateBatch (batchSize*getBatchFraction());
   }

   /**
    * Simulates a batch with length `batchLength`. This method initializes
    * the model-specific counters by using  #initBatchStat, simulates the
    * batch, and adds observations using  #addRealBatchObs. It also calls
    * #addEffectiveBatchObs if aggregation is turned OFF. If the batch
    * lengths are not fixed, one can call this method to call
    * #simulateBatch(double) with <tt>Double.POSITIVE_INFINITY</tt>; this
    * prevents the method from scheduling the ending event, and let the
    * simulator call `simulator().stop` at appropriate time.
    */
   public void simulateBatch (double batchLength) {
      initBatchStat();
      endSimEvent.schedule (batchLength);
      simulator().start();
      endSimEvent.cancel();
      batchFraction = 1.0;
      batchSizeMultiplier = 1.0;
      aggregateUpdated = false;
      ++doneBatches;
      if (batchLengths || aggregation)
         batchTimes.add (simulator().time());
      else {
         batchTimes.setQuick (0, batchTimes.getQuick (1));
         batchTimes.setQuick (1, simulator().time());
      }
      addRealBatchObs();
      if (!aggregation) {
         nAgr = 1;
         addEffectiveBatchObs (doneBatches - 1, 1, getRealBatchLength (doneBatches - 1));
      }
   }

   /**
    * Adjusts the target number of real batches to simulate
    * `numNewBatches` additionnal real batches. This method clamps the
    * target number of real batches to the maximal number of batches, and
    * ensures that the target number of batches is a multiple of the
    * minimal number of batches when aggregation is turned ON.
    */
   public void adjustTargetBatches (int numNewBatches) {
      if (doneBatches + numNewBatches > maxBatches) {
         numNewBatches = maxBatches - doneBatches;
         if (aggregation && numNewBatches < minBatches)
            // We must simulate at least minBatches real batches, so
            // do nothing
            return;
      }
      if (aggregation) {
         // For effective batches to change,
         // we need to simulate at least minBatches additional real batches.
         // The number of new batches is therefore rounded to
         // the smallest multiple greater than or
         // equal to minBatches.
         if (numNewBatches % minBatches != 0) {
            numNewBatches /= minBatches;
            numNewBatches *= minBatches;
            if (doneBatches + numNewBatches + minBatches <= maxBatches)
               numNewBatches += minBatches;
         }
      }
      targetBatches = doneBatches + numNewBatches;
   }

   /**
    * Simulates batches until the number of completed real batches
    * corresponds to the target number of batches. This method first
    * allocates the capacity for simulating  #getTargetBatches. It then
    * simulate each batch using  #simulateBatch. If aggregation is turned
    * ON, this method also calls  #initEffectiveBatchProbes, and
    * #addEffectiveBatchObs to manage effective batches.
    */
   public void simulateBatches () {
      if (doneBatches < targetBatches)
         // Increased the target number of required batches
         // and allocate memory to store the observations
         // if this is required.
         allocateCapacity();
      while (doneBatches < targetBatches)
         // targetBatches may be changed by a simulation event.
         // This can be used to stop the simulation before the
         // desired number of batches is obtained, e.g., when
         // the simulated system seems unstable.
         simulateBatch();
      if (aggregation)
         makeAggregateBatches();
   }

   /**
    * Performs a batch means simulation. This method resets the state of
    * the system by calling  #init, and calls  #warmup to perform the
    * warmup. Then, the method calls  #simulateBatches and
    * #getRequiredNewBatches until the number of completed real batches
    * equals or exceeds the target number of batches.
    */
   public void simulate() {
      init();
      simulating = true;
      try {
         warmup (warmupTime);
         while (doneBatches < targetBatches) {
            simulateBatches();
            adjustTargetBatches (getRequiredNewBatches());
         }
         if (aggregation && !aggregateUpdated)
            makeAggregateBatches();
      }
      finally {
         simulating = false;
      }
   }


   private final void makeAggregateBatches() {
      // When using fixed number of batches, we must
      // add the observations to statistical collectors
      // only at the end of the simulation.
      initEffectiveBatchProbes();
      // nAgr gives the number of real batches in each effective batch
      nAgr = (doneBatches - droppedBatches) / minBatches;
      // tnb gives the number of observations in the statistical counters
      int tnb = minBatches;
      if (nAgr == 0) {
         // doneBatches is smaller than minBatches, but we would like
         // observations to be collected and displayed.
         // We will have less than minBatches observations.
         nAgr = 1;
         tnb = doneBatches - droppedBatches;
      }
      for (int i = 0; i < tnb; i++) {
         int startAgr = i*nAgr + droppedBatches;
         addEffectiveBatchObs
            (startAgr, nAgr, batchTimes.getQuick (startAgr - droppedBatches + nAgr)
             - batchTimes.getQuick (startAgr - droppedBatches));
      }
      aggregateUpdated = true;
   }

   private static final class EndSimEvent extends Event {
      public EndSimEvent (Simulator sim) {
         super (sim);
      }

      public void actions() {
         simulator().stop();
      }
   }

   public String toString() {
      StringBuffer sb = new StringBuffer (getClass().getName());
      sb.append ('[');
      sb.append ("minimal number of batches: ").append (minBatches);
      if (maxBatches < Integer.MAX_VALUE)
         sb.append (", maximal number of batches: ").append (maxBatches);
      sb.append (", target number of batches: ").append (targetBatches);
      sb.append (", warmup time: ").append (warmupTime);
      if (simulating)
         sb.append (", simulation in progress");
      else
         sb.append (", simulation stopped");
      if (warmupDone)
         sb.append (", number of completed batches: ").append (doneBatches);
      else
         sb.append (", warmup not completed");
      if (aggregation)
         sb.append (", batch aggregation turned ON");
      else {
         sb.append (", batch aggregation turned OFF");
         if (batchLengths)
            sb.append (", batch lengths are kept");
         else
            sb.append (", batch lengths are not kept");
      }
      sb.append (']');
      return sb.toString();
   }

/**
 * Returns the sum of elements `start`, …, `start + length - 1`, in the array
 * `a`.
 *  @param a            the source array.
 *  @param start        the index of the first element to sum.
 *  @param length       the number of elements in the sum.
 *  @return the sum of the elements.
 *
 *  @exception NullPointerException if `a` is `null`.
 *  @exception IllegalArgumentException `length` is negative.
 *  @exception IndexOutOfBoundsException if `start` is negative, or `start +
 * length` is greater than the length of `a`.
 */
public static double getSum (double[] a, int start, int length) {
      if (start < 0)
         throw new IndexOutOfBoundsException ("start must not be negative");
      if (length < 0)
         throw new IllegalArgumentException ("length is negative");
      if (start + length > a.length)
         throw new IndexOutOfBoundsException ("Not enough elements in the array");
      double s = 0;
      for (int i = start; i < start + length; i++)
         s += a[i];
      return s;
   }

   /**
    * Returns the sum of elements `start`, …, `start + length - 1`, in the
    * array list `l`.
    *  @param l            the source array list.
    *  @param start        the index of the first element to sum.
    *  @param length       the number of elements in the sum.
    *  @return the sum of the elements.
    *
    *  @exception NullPointerException if `l` is `null`.
    *  @exception IllegalArgumentException `length` is negative.
    *  @exception IndexOutOfBoundsException if `start` is negative, or
    * `start + length` is greater than the size of `l`.
    */
   public static double getSum (DoubleArrayList l, int start, int length) {
      if (start < 0)
         throw new IndexOutOfBoundsException ("start must not be negative");
      if (length < 0)
         throw new IllegalArgumentException ("length is negative");
      if (start + length > l.size())
         throw new IndexOutOfBoundsException ("Not enough elements in the array list");
      double s = 0;
      for (int i = start; i < start + length; i++)
         s += l.getQuick (i);
      return s;
   }

   /**
    * Returns the sum of elements `start`, …, `start + length - 1`, in the
    * 1D matrix `m`.
    *  @param m            the source 1D matrix.
    *  @param start        the index of the first element to sum.
    *  @param length       the number of elements in the sum.
    *  @return the sum of the elements.
    *
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException `length` is negative.
    *  @exception IndexOutOfBoundsException if `start` is negative, or
    * `start + length` is greater than the size of `m`.
    */
   public static double getSum (DoubleMatrix1D m, int start, int length) {
      if (start < 0)
         throw new IndexOutOfBoundsException ("start must not be negative");
      if (length < 0)
         throw new IllegalArgumentException ("length is negative");
      if (start + length > m.size())
         throw new IndexOutOfBoundsException ("Not enough elements in the matrix");
      double s = 0;
      for (int i = start; i < start + length; i++)
         s += m.getQuick (i);
      return s;
   }

   /**
    * Returns an array containing the sum of columns `startColumn`, …,
    * `startColumn + numColumns - 1`, in the 2D matrix represented by the
    * 2D array `a`. The given array is assumed to be rectangular, i.e.,
    * each of its array elements has the same length.
    *  @param a            the source 2D array.
    *  @param startColumn  the index of the first column to sum.
    *  @param numColumns   the number of columns in the sum.
    *  @return an array of `a.length` elements containing the sums of the
    * columns.
    *
    *  @exception NullPointerException if `a` is `null`.
    *  @exception IllegalArgumentException `numColumns` is negative.
    *  @exception IndexOutOfBoundsException if `startColumn` is negative,
    * or `startColumn + numColumns` is greater than the number of columns
    * in the matrix represented by `a`.
    */
   public static double[] getSum (double[][] a, int startColumn,
                                  int numColumns) {
      if (a.length == 0)
         return new double[0];
      if (startColumn < 0)
         throw new IndexOutOfBoundsException ("startColumn must not be negative");
      if (numColumns < 0)
         throw new IllegalArgumentException ("numColumns is negative");
      if (startColumn + numColumns > a[0].length)
         throw new IndexOutOfBoundsException ("Not enough elements in the array");
      double[] res = new double[a.length];
      for (int i = 0; i < res.length; i++)
         for (int j = startColumn; j < startColumn + numColumns; j++)
            res[i] += a[i][j];
      return res;
   }

   /**
    * Returns an array containing the sum of columns `startColumn`, …,
    * `startColumn + numColumns - 1`, in the 2D matrix `m`.
    *  @param m            the source 2D matrix.
    *  @param startColumn  the index of the first column to sum.
    *  @param numColumns   the number of columns in the sum.
    *  @return an array of `m.rows()` elements containing the sum of the
    * columns.
    *
    *  @exception NullPointerException if `m` is `null`.
    *  @exception IllegalArgumentException `numColumns` is negative.
    *  @exception IndexOutOfBoundsException if `startColumn` is negative,
    * or `startColumn + numColumns` is greater than the number of columns
    * in `m`.
    */
   public static double[] getSum (DoubleMatrix2D m, int startColumn,
                                  int numColumns) {
      if (startColumn < 0)
         throw new IndexOutOfBoundsException ("startColumn must not be negative");
      if (numColumns < 0)
         throw new IllegalArgumentException ("numColumns is negative");
      if (startColumn + numColumns > m.columns())
         throw new IndexOutOfBoundsException ("Not enough columns in the matrix");
      double[] res = new double[m.rows()];
      for (int i = 0; i < res.length; i++)
         for (int j = startColumn; j < startColumn + numColumns; j++)
            res[i] += m.getQuick (i, j);
      return res;
   }

   /**
    * Regroups the elements in array `a` by summing each successive `x`
    * values. When this method returns, element `i` of the given array
    * corresponds to the sum of elements `ix`, …, `ix + x - 1` in the
    * original array. If the size of the array is not a multiple of `x`,
    * the remaining elements are summed up and added into an extra element
    * of the transformed array. Remaining elements of the transformed
    * array are set to 0.
    *  @param a            the array being processed.
    *  @param x            the regrouping factor.
    *  @exception NullPointerException if `a` is `null`.
    *  @exception IllegalArgumentException if `x` is smaller than 1.
    */
   public static void regroupElements (double[] a, int x) {
      if (x < 1)
         throw new IllegalArgumentException ("x < 1");
      int gs = a.length / x;
      for (int i = 0; i < gs; i++) {
         double o = 0;
         for (int j = 0; j < x; j++)
            o += a[i*x + j];
         a[i] = o;
      }
      int m = a.length % x;
      if (m > 0) {
         double o = 0;
         int r = a.length - x*gs;
         for (int j = 0; j < r; j++)
            o += a[gs*x + j];
         a[gs] = o;
         ++gs;
      }
      for (int i = gs; i < a.length; i++)
         a[i] = 0;
   }

   /**
    * Same as  #regroupElements(double[],int) for an array list. The size
    * of the list is also divided by `x`.
    *  @param l            the array list being processed.
    *  @param x            the regrouping factor.
    *  @exception NullPointerException if `l` is `null`.
    *  @exception IllegalArgumentException if `x` is smaller than 1.
    */
   public static void regroupElements (DoubleArrayList l, int x) {
      if (x < 1)
         throw new IllegalArgumentException ("x < 1");
      int gs = l.size() / x;
      for (int i = 0; i < gs; i++) {
         double o = 0;
         for (int j = 0; j < x; j++)
            o += l.getQuick (i*x + j);
         l.setQuick (i, o);
      }
      int m = l.size() % x;
      if (m > 0) {
         double o = 0;
         int r = l.size() - x*gs;
         for (int j = 0; j < r; j++)
            o += l.getQuick (gs*x + j);
         l.setQuick (gs, o);
         ++gs;
      }
      l.setSize (gs);
   }

   /**
    * Same as  #regroupElements(double[],int) for a 1D matrix.
    *  @param mat          the 1D matrix being processed.
    *  @param x            the regrouping factor.
    *  @exception NullPointerException if `mat` is `null`.
    *  @exception IllegalArgumentException if `x` is smaller than 1.
    */
   public static void regroupElements (DoubleMatrix1D mat, int x) {
      if (x < 1)
         throw new IllegalArgumentException ("x < 1");
      int gs = mat.size() / x;
      for (int i = 0; i < gs; i++) {
         double o = 0;
         for (int j = 0; j < x; j++)
            o += mat.getQuick (i*x + j);
         mat.setQuick (i, o);
      }
      int m = mat.size() % x;
      if (m > 0) {
         double o = 0;
         int r = mat.size() - x*gs;
         for (int j = 0; j < r; j++)
            o += mat.getQuick (gs*x + j);
         mat.setQuick (gs, o);
         ++gs;
      }
      for (int i = gs; i < mat.size(); i++)
         mat.setQuick (i, 0);
   }

   /**
    * Same as  #regroupElements(double[],int) for a 2D matrix. This method
    * regroups columns and considers each row as an independent array.
    *  @param mat          the 2D matrix being processed.
    *  @param x            the regrouping factor.
    *  @exception NullPointerException if `mat` is `null`.
    *  @exception IllegalArgumentException if `x` is smaller than 1.
    */
   public static void regroupElements (DoubleMatrix2D mat, int x) {
      if (x < 1)
         throw new IllegalArgumentException ("x < 1");
      int gs = mat.columns() / x;
      for (int i = 0; i < gs; i++) {
         for (int r = 0; r < mat.rows(); r++) {
            double o = 0;
            for (int j = 0; j < x; j++)
               o += mat.getQuick (r, i*x + j);
            mat.setQuick (r, i, o);
         }
      }
      int m = mat.columns() % x;
      if (m > 0) {
         int r = mat.columns() - x*gs;
         for (int row = 0; row < mat.rows(); row++) {
            double o = 0;
            for (int j = 0; j < r; j++)
               o += mat.getQuick (row, gs*x + j);
            mat.setQuick (row, gs, o);
         }
         ++gs;
      }
      for (int row = 0; row < mat.rows(); row++)
         for (int i = gs; i < mat.columns(); i++)
            mat.setQuick (row, i, 0);
   }
}