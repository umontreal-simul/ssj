/*
 * Class:        VarianceGammaProcess
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since        2004
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
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;

/**
 * This class represents a *variance gamma* (VG) process @f$\{S(t) = X(t;
 * \theta, \sigma, \nu) : t \geq0\}@f$. This process is obtained as a
 * subordinate of the Brownian motion process @f$B(t;\theta,\sigma)@f$
 * using the operational time @f$G(t;1,\nu)@f$ (see @cite pFEL66a,
 * @cite fAVR06a&thinsp;):
 * @anchor REF_stochprocess_VarianceGammaProcess_VGeqn
 * @f[
 *   X(t; \theta, \sigma, \nu) := B(G(t;1,\nu),\theta, \sigma). \tag{VGeqn}
 * @f]
 * See also @cite fMAD98a, @cite fMAD91a, @cite fMAD90a&thinsp; for
 * applications to modelling asset returns and option pricing.
 *
 * The process is sampled as follows: when `generatePath()` is called, the
 * method <tt>generatePath()</tt> of the inner  @ref GammaProcess is called;
 * its path is then used to set the observation times of the
 * @ref BrownianMotion. Finally, the method `generatePath()` of the
 * @ref BrownianMotion is called. *Warning*: If one wants to reduced the
 * variance as much as possible in a QMC simulation, this way of proceeding
 * is not optimal. Use the method `generatePath(uniform01)` instead.
 *
 * If one calls the `nextObservation` method, the operational time is
 * generated first, followed by the corresponding brownian motion increment,
 * which is then returned.
 *
 * Note that if one wishes to use *bridge* sampling with the
 * `nextObservation` method, both the gamma process @f$G@f$ and the Brownian
 * motion process @f$B@f$ should use bridge sampling so that their
 * observations are synchronized.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class VarianceGammaProcess extends StochasticProcess {

    protected GammaProcess   randomTime;  // For the transformed time method
    protected BrownianMotion BM;

    protected double       theta,
                           sigma,
                           nu;


   public VarianceGammaProcess() {}

/**
 * Constructs a new `VarianceGammaProcess` with parameters @f$\theta=
 * \mathtt{theta}@f$, @f$\sigma= \mathtt{sigma}@f$, @f$\nu= \mathtt{nu}@f$
 * and initial value @f$S(t_0) = \mathtt{s0}@f$. `stream` is used to generate
 * both the  @ref BrownianMotion @f$B@f$ and the  @ref GammaProcess @f$G@f$
 * in ( {@link REF_stochprocess_VarianceGammaProcess_VGeqn
 * VGeqn} ).
 */
public VarianceGammaProcess (double s0, double theta, double sigma,
                                double nu, RandomStream stream) {
        this (s0, new BrownianMotion (s0, theta, sigma, stream),
                  new GammaProcess (0.0, 1.0, nu, stream));
    }

   /**
    * Constructs a new `VarianceGammaProcess`. The parameters
    * @f$\theta@f$ and @f$\sigma@f$ are set to the parameters @f$\mu@f$
    * and @f$\sigma@f$, respectively, of the  @ref BrownianMotion `BM`
    * and the parameter @f$\nu@f$ is set to the parameter @f$\nu@f$ of
    * the  @ref GammaProcess `Gamma`. The parameters @f$\mu@f$ and
    * @f$x0@f$ of the  @ref GammaProcess are overwritten to equal 1 and 0
    * respectively. The initial value of the process is @f$S(t_0) =
    * {\mathtt{s0}}@f$.
    */
   public VarianceGammaProcess (double s0, BrownianMotion BM,
                                GammaProcess Gamma) {
        this.BM         = BM;
        Gamma.setParams(0.0, 1.0, Gamma.getNu());  // forces the average of the GammaProcess
        randomTime      = Gamma;                   // to be 1.0 and the initial value to be 0.0
        setParams (s0, BM.getMu(), BM.getSigma(), Gamma.getNu());
    }

   /**
    * Generates the observation for the next time. It also works with
    * bridge sampling; however *both*  @ref BrownianMotionBridge and
    * @ref GammaProcessBridge must be used in the constructor in that
    * case. Furthermore, for bridge sampling, the order of the
    * observations is that of the bridge, not sequential order.
    */
   public double nextObservation() {
        // We first generate w, then verify what its new counter value is
        // This is necessary to be general enough to handle bridge sampling
        double nextBM = BM.nextObservation (randomTime.nextObservation ());
        observationIndex = BM.getCurrentObservationIndex();
        path[observationIndex] = nextBM;
        observationCounter++;
        return nextBM;
    }

   /**
    * Generates and returns the path. To do so, it first generates the
    * complete path of the inner  @ref GammaProcess and sets the
    * observation times of the inner  @ref BrownianMotion to this path.
    * This method is not optimal to reduce the variance in QMC
    * simulations; use `generatePath(double[] uniform01)` for that.
    */
   public double[] generatePath() {
        BM.setObservationTimes(randomTime.generatePath(), d);
        path = BM.generatePath();
        observationIndex = d;
        observationCounter = d;
        return path;
    }

   /**
    * Similar to the usual `generatePath()`, but here the uniform random
    * numbers used for the simulation must be provided to the method. This
    * allows to properly use the uniform random variates in QMC
    * simulations. This method divides the table of uniform random numbers
    * `uniform01` in two smaller tables, the first one, containing the odd
    * indices of `uniform01` which are used to generate the path of the
    * inner  @ref GammaProcess, and the even indices (in the second table)
    * are used to generate the path of the inner  @ref BrownianMotion.
    * This way of proceeding reduces the variance as much as possible for
    * QMC simulations.
    */
   public double[] generatePath (double[] uniform01) {
        int dd = uniform01.length;
        int d = dd / 2;

        if (dd % 2 != 0) {
            throw new IllegalArgumentException (
                     "The Array uniform01 must have a even length");
        }

        double[] QMCpointsGP = new double[d];
        double[] QMCpointsBM = new double[d];

        for(int i = 0; i < d; i++){
            QMCpointsGP[i] = uniform01[2 * i];  // the odd numbers for the gamma process
            QMCpointsBM[i] = uniform01[2 * i + 1];  //  and the even for the BM process
        }
        BM.setObservationTimes(randomTime.generatePath(QMCpointsGP), d);

        path = BM.generatePath(QMCpointsBM);
        observationIndex = d;
        observationCounter = d;
        return path;
    }

    /**
     * Resets the observation index and counter to 0 and applies the
     * `resetStartProcess` method to the  @ref BrownianMotion and the
     * @ref GammaProcess objects used to generate this process.
     */
    public void resetStartProcess() {
        observationIndex   = 0;
        observationCounter = 0;
        BM.resetStartProcess();
        randomTime.resetStartProcess();
    }

   /**
    * Sets the parameters @f$S(t_0) =@f$ `s0`, @f$\theta=@f$ `theta`,
    * @f$\sigma=@f$ `sigma` and @f$\nu=@f$ `nu` of the process.
    * *Warning*: This method will recompute some quantities stored
    * internally, which may be slow if called repeatedly.
    */
   public void setParams (double s0, double theta, double sigma, double nu) {
        this.x0    = s0;
        this.theta = theta;
        this.sigma = sigma;
        this.nu    = nu;
        if (observationTimesSet) init(); // Otherwise no need to.
    }

   /**
    * Returns the value of the parameter @f$\theta@f$.
    */
   public double getTheta() { return theta; }

   /**
    * Returns the value of the parameter @f$\sigma@f$.
    */
   public double getSigma() { return sigma; }

   /**
    * Returns the value of the parameter @f$\nu@f$.
    */
   public double getNu() { return nu; }


    protected void init() {
        super.init();
        if(observationTimesSet){
            randomTime.setObservationTimes(t, d);
            randomTime.x0 = t[0];
        }
    }

/**
 * Sets the observation times on the `VarianceGammaProcess` as usual, but
 * also sets the observation times of the underlying  @ref GammaProcess. It
 * furthermore sets the starting *value* of the  @ref GammaProcess to `t[0]`.
 */
public void setObservationTimes (double t[], int d) {
        super.setObservationTimes(t, d);  //sets the observation times of the GammaProcess by
    }                                     //calling init()

    /**
     * Resets the  @ref umontreal.ssj.rng.RandomStream ’s. Warning: this
     * method sets both the  @ref umontreal.ssj.rng.RandomStream of the
     * @ref BrownianMotion and of the  @ref GammaProcess to the same
     * @ref umontreal.ssj.rng.RandomStream.
     */
    public void setStream (RandomStream stream) {
         BM.setStream (stream);
   }

   /**
    * Returns the random stream of the  @ref BrownianMotion process, which
    * should be the same as for the  @ref GammaProcess.
    */
   public RandomStream getStream() {
        return BM.getStream();
   }

   /**
    * Returns a reference to the inner  @ref BrownianMotion.
    */
   public BrownianMotion getBrownianMotion() {
      return BM;
   }

   /**
    * Returns a reference to the inner  @ref GammaProcess.
    */
   public GammaProcess getGammaProcess() {
      return randomTime;
   }

}