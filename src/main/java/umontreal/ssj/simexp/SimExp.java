/*
 * Class:        SimExp
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

import umontreal.ssj.simevents.Simulator;
import umontreal.ssj.stat.FunctionOfMultipleMeansTally;
import umontreal.ssj.stat.StatProbe;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stat.matrix.MatrixOfStatProbes;

/**
 * Represents a framework for performing experiments using simulation. This
 * class defines an abstract  #simulate method that should implement the
 * simulation logic. It also provides utility methods to estimate the
 * required number of additional observations that would be necessary for an
 * estimator to reach a given precision, for sequential sampling.
 *
 * This class is the base class of  @ref BatchMeansSim and  @ref RepSim
 * implementing the logic for a simulation on infinite and finite horizon,
 * respectively.
 *
 * <div class="SSJ-bigskip"></div>
 */
public abstract class SimExp {
   protected Simulator sim;

/**
 * Determines if the simulation is in progress.
 */
protected boolean simulating = false;

   /**
    * Constructs a new object for performing experiments using the default
    * simulator returned by `Simulator.getDefaultSimulator()`.
    */
   protected SimExp () {
      this (Simulator.getDefaultSimulator());
   }

   /**
    * Constructs a new object performing experiments using the given
    * simulator `sim`.
    *  @param sim          the simulator attached to this object.
    */
   protected SimExp (Simulator sim) {
      if (sim == null)
         throw new NullPointerException();
      this.sim = sim;
   }

   /**
    * Returns the simulator linked to this experiment object.
    *  @return the simulator linked to the experiment object
    */
   public final Simulator simulator () {
      return sim;
   }

   /**
    * Sets the simulator associated with this experiment to `sim`. This
    * method should not be called while this object is simulating.
    *  @param sim          the new simulator
    */
   public final void setSimulator (Simulator sim) {
       if (sim == null)
          throw new NullPointerException();
      if (simulating)
         throw new UnsupportedOperationException("Unable to set Simulator, experiment object already running");
      this.sim = sim;
   }

   /**
    * Determines if the simulation is in progress.
    *  @return `true` if and only if simulation is in progress.
    */
   public boolean isSimulating() {
      return simulating;
   }

   /**
    * Performs an experiment whose logic depends on the used subclass.
    * Before starting the simulation, this method should set  #simulating
    * to `true`, and reset it to `false` after the simulation is done. It
    * is recommended to reset  #simulating to `false` inside a `finally`
    * clause to prevent the indicator from remaining `true` in the case of
    * error during simulation.
    *  @exception IllegalStateException if  #simulating is `true` when
    * calling this method.
    */
   public abstract void simulate();

   /**
    * Returns the approximate number of additional observations required
    * to reach a relative error smaller than or equal to `targetError` for
    * each tally in the array `a` when confidence intervals are computed
    * with confidence level `level`. For each statistical collector in the
    * given array, a confidence interval is computed independently of the
    * other collectors, and an error check is performed by
    * #getRequiredNewObservations(StatProbe,double,double) to determine
    * the required number of additional observations. The method returns
    * the maximal number of required observations.
    *  @param a            the array of probes.
    *  @param targetError  the target relative error.
    *  @param level        the desired probability that, for a given
    *                      statistical collector, the (random) confidence
    *                      interval covers the true mean (a constant).
    *  @return an estimate of the required number of additional
    * observations to reach the precision.
    */
   public static int getRequiredNewObservations (StatProbe[] a,
                                                 double targetError,
                                                 double level) {
      int nnewobs = 0;
      for (StatProbe sp : a) {
         int re = sp == null ? 0 : getRequiredNewObservations (sp, targetError, level);
         if (re > nnewobs)
            nnewobs = re;
      }
      return nnewobs;
   }

   /**
    * Returns the approximate number of additional observations required
    * to reach a relative error smaller than or equal to `targetError` for
    * each tally enumerated by `it` when confidence intervals are computed
    * with confidence level `level`. For each statistical collector
    * returned by the iterator obtained from `it`, a confidence interval
    * is computed independently of the other collectors, and an error
    * check is performed by
    * #getRequiredNewObservations(StatProbe,double,double) to determine
    * the required number of additional observations. The method returns
    * the maximal number of required observations.
    *  @param it           the iterable used to enumerate probes.
    *  @param targetError  the target relative error.
    *  @param level        the desired probability that, for a given
    *                      statistical collector, the (random) confidence
    *                      interval covers the true mean (a constant).
    *  @return an estimate of the required number of additional
    * observations to reach the precision.
    */
   public static int getRequiredNewObservations (
                                     Iterable<? extends StatProbe> it,
                                                 double targetError,
                                                 double level) {
      int nnewobs = 0;
      for (StatProbe sp : it) {
         int re = sp == null ? 0 : getRequiredNewObservations (sp, targetError, level);
         if (re > nnewobs)
            nnewobs = re;
      }
      return nnewobs;
   }

   /**
    * Calls  #getRequiredNewObservations(double,double,int,double) with
    * the average, confidence interval radius, and number of observations
    * given by the statistical probe `probe`. This method always returns 0
    * if the probe is not a tally. For a  @ref umontreal.ssj.stat.Tally,
    * the confidence interval is computed using
    * umontreal.ssj.stat.Tally.confidenceIntervalStudent(double,double[]).
    * For a  @ref umontreal.ssj.stat.FunctionOfMultipleMeansTally, it is
    * computed using
    * umontreal.ssj.stat.FunctionOfMultipleMeansTally.confidenceIntervalDelta(double,double[]).
    *  @param probe        the statistical probe being checked.
    *  @param targetError  the target relative error.
    *  @param level        the desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant).
    *  @return the number of required additional observations.
    */
   public static int getRequiredNewObservations (StatProbe probe,
                                                 double targetError,
                                                 double level) {
      if (probe instanceof Tally)
         return getRequiredNewObservationsTally ((Tally) probe, targetError, level);
      else if (probe instanceof FunctionOfMultipleMeansTally)
         return getRequiredNewObservationsTally ((FunctionOfMultipleMeansTally) probe, targetError, level);
      else
         return 0;
   }

   /**
    * Calls  #getRequiredNewObservations(double,double,int,double) with
    * the average, confidence interval radius, and number of observations
    * given by the tally `ta`. The confidence interval is computed using
    * umontreal.ssj.stat.Tally.confidenceIntervalStudent(double,double[]).
    *  @param ta           the tally being checked.
    *  @param targetError  the target relative error.
    *  @param level        the desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant).
    *  @return the number of required additional observations.
    */
   public static int getRequiredNewObservationsTally (Tally ta,
                                                 double targetError,
                                                 double level) {
      double[] cr = new double[2];
      int no = ta.numberObs();
      if (no >= 2)
         ta.confidenceIntervalStudent (level, cr);

      return getRequiredNewObservations (cr[0], cr[1], no, targetError);
   }

   /**
    * Calls  #getRequiredNewObservations(double,double,int,double) with
    * the average, confidence interval radius, and number of observations
    * given by the function of multiple means `fmmt`. The confidence
    * interval is computed using
    * umontreal.ssj.stat.FunctionOfMultipleMeansTally.confidenceIntervalDelta(double,double[]).
    *  @param fmmt         the function of multiple means being checked.
    *  @param targetError  the target relative error.
    *  @param level        the desired probability that the (random)
    *                      confidence interval covers the true mean (a
    *                      constant).
    *  @return the number of required additional observations.
    */
   public static int getRequiredNewObservationsTally (
                                      FunctionOfMultipleMeansTally fmmt,
                                                 double targetError,
                                                 double level) {
      double[] cr = new double[2];
      int no = fmmt.numberObs();
      if (no >= 2)
         fmmt.confidenceIntervalDelta (level, cr);

      return getRequiredNewObservations (cr[0], cr[1], no, targetError);
   }

   /**
    * Returns the approximate number of additional observations needed for
    * the point estimator @f$\bar{X}_n=@f$&nbsp;`center`, computed using
    * @f$n=@f$&nbsp;`numberObs` observations and with a confidence
    * interval having radius @f$\delta_n/\sqrt{n}=@f$&nbsp;`radius`, to
    * have a relative error less than or equal to
    * @f$\epsilon=@f$&nbsp;`targetError`. It is assumed that
    * @f$\bar{X}_n@f$ is an estimator of a mean @f$\mu@f$, @f$n@f$ is the
    * number of observations `numberObs`, and that
    * @f$\delta_n/\sqrt{n}\to0@f$ when @f$n\to\infty@f$.
    *
    * If @f$n@f$ is less than 1, this method returns 0. Otherwise, the
    * relative error given by @f$\delta_n/|\sqrt{n}\bar{X}_n|@f$ should
    * be smaller than or equal to @f$\epsilon@f$. If the inequality is
    * true, this returns 0. Otherwise, the minimal @f$n^*@f$ for which
    * this inequality holds is approximated as follows. The target radius
    * is given by @f$\delta^*=\epsilon|\mu|@f$, which is approximated
    * by @f$\epsilon|\bar{X}_n|<\delta_n/\sqrt{n}@f$. The method must
    * select @f$n^*@f$ for which
    * @f$\delta_{n^*}/\sqrt{n^*}\le\delta^*@f$, which will be
    * approximately true if
    * @f$\delta_{n^*}/\sqrt{n^*}\le\epsilon|\bar{X}_n|@f$. Therefore,
    * @f[
    *   n^*\ge(\delta_{n^*}/(\epsilon|\bar{X}_n|))^2\approx(\delta_n/(\epsilon|\bar{X}_n|))^2.
    * @f]
    * The method returns
    * @f$\mathrm{round}((\delta_n\sqrt{n}/(\epsilon|\bar{X}_n|))^2)-n@f$
    * where @f$\mathrm{round}(\cdot)@f$ rounds its argument to the
    * nearest integer.
    *  @param center       the value of the point estimator.
    *  @param radius       the radius of the confidence interval.
    *  @param numberObs    the number of observations.
    *  @param targetError  the target relative error.
    *  @return an estimate of the required number of additional
    * observations to reach the precision.
    *
    *  @exception IllegalArgumentException if `radius` or `targetError`
    * are negative.
    */
   public static int getRequiredNewObservations (double center,
                                                 double radius,
                                                 int numberObs,
                                                 double targetError) {
      if (radius < 0)
         throw new IllegalArgumentException
            ("The radius must not be negative");
      if (targetError < 0)
         throw new IllegalArgumentException
            ("The target error must not be negative");
      if (numberObs < 1)
         return 0;
      double targetRadius = targetError*Math.abs (center);
      if (radius <= targetRadius)
         return 0;
      double deltan = radius*Math.sqrt (numberObs);
      double sqrtnp = deltan/targetRadius;
      double totalnobs = sqrtnp*sqrtnp;
      int nn = (int)Math.round (totalnobs) - numberObs;
      if (nn < 0)
         return 0;
      return nn;
   }

}