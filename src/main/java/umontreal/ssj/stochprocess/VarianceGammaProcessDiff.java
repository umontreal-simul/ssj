/*
 * Class:        VarianceGammaProcessDiff
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
 * \theta, \sigma, \nu) : t \geq0\}@f$. This process is generated using
 * *difference of gamma sampling* (see @cite fAVR03a, @cite fAVR06a&thinsp;),
 * which uses the representation of the VG process as the difference of two
 * independent  @ref GammaProcess ’es (see @cite fMAD98a&thinsp;):
 * @anchor REF_stochprocess_VarianceGammaProcessDiff_dblGammaEqn
 * @f[
 *   X(t; \theta, \sigma, \nu) := X(0) + \Gamma^+(t; \mu_p, \nu_p) - \Gamma^-(t; \mu_n, \nu_n) \tag{dblGammaEqn}
 * @f]
 * where @f$X(0)@f$ is a constant corresponding to the initial value of the
 * process and
 * @anchor REF_stochprocess_VarianceGammaProcessDiff_dblGammaParams
 * @f[
 *   \begin{array}{rcl}
 *    \mu_p 
 *    & 
 *    = 
 *    & 
 *    (\sqrt{ \theta^2 + 2\sigma^2/\nu} + \theta)/2 
 *    \\ 
 *    \mu_n 
 *    & 
 *    = 
 *    & 
 *    (\sqrt{ \theta^2 + 2\sigma^2/\nu} - \theta)/2 
 *    \\ 
 *    \nu_p 
 *    & 
 *    = 
 *    & 
 *    \nu\mu_p^2
 *    \\ 
 *    \nu_n 
 *    & 
 *    = 
 *    & 
 *    \nu\mu_n^2 
 *   \end{array} \tag{dblGammaParams}
 * @f]
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class VarianceGammaProcessDiff extends VarianceGammaProcess {
    protected GammaProcess gpos;
    protected GammaProcess gneg;
    protected double       mup, mun,
                           nup, nun;

   /**
    * Constructs a new `VarianceGammaProcessDiff` with parameters
    * @f$\theta= \mathtt{theta}@f$, @f$\sigma= \mathtt{sigma}@f$,
    * @f$\nu= \mathtt{nu}@f$ and initial value @f$S(t_0) =
    * \mathtt{s0}@f$. `stream` is used by two instances of
    * @ref GammaProcess, @f$\Gamma^+@f$ and @f$\Gamma^-@f$,
    * respectively. The other parameters are as in the class
    * @ref VarianceGammaProcess. The  @ref GammaProcess objects for
    * @f$\Gamma^+@f$ and @f$\Gamma^-@f$ are constructed using the
    * parameters from (
    * {@link REF_stochprocess_VarianceGammaProcessDiff_dblGammaParams
    * dblGammaParams} ) and their initial values @f$\Gamma^+(t_0)@f$ and
    * @f$\Gamma^-(t_0)@f$ are set to @f$0@f$.
    */
   public VarianceGammaProcessDiff (double s0, double theta, double sigma,
                                    double nu, RandomStream stream) {
        this (s0, theta, sigma, nu,
              new GammaProcess (0.0, 1.0, 1.0, stream),
              new GammaProcess (0.0, 1.0, 1.0, stream));
        // Params mu, nu of the 2 gamma processes are redefined in init()
        // which will be called after a call to 'setObservTimes'
    }

   /**
    * The parameters of the  @ref GammaProcess objects for @f$\Gamma^+@f$
    * and @f$\Gamma^-@f$ are set to those of (
    * {@link REF_stochprocess_VarianceGammaProcessDiff_dblGammaParams
    * dblGammaParams} ) and their initial values @f$\Gamma^+(t_0)@f$ and
    * @f$\Gamma^-(t_0)@f$ are set to @f$t_0@f$. The `RandomStream` of the
    * @f$\Gamma^-@f$ process is overwritten with the `RandomStream` of
    * the @f$\Gamma^+@f$ process.
    */
   public VarianceGammaProcessDiff (double s0, double theta, double sigma,
                                    double nu, GammaProcess gpos,
                                    GammaProcess gneg) {
        this.gpos = gpos;
        this.gneg = gneg;
        setParams (s0, theta, sigma, nu);
        gneg.setStream(gpos.getStream());  // to avoid confusion with stream because there is only
                                 // one stream in the other constructor
    }


   public double nextObservation() {
        // This implementation takes possible bridge sampling into account
        double s = x0 + gpos.nextObservation() - gneg.nextObservation();
        observationIndex = gpos.getCurrentObservationIndex();
        path[observationIndex] = s;
        observationCounter++;
        return s;
     }

// no longer useful, this method was created to automaticaly alternate
// between the two processes the uniform random variables used in the
// in the simulation.  However, this method does not work if the two
// GammaProcess are PCA...
//    public double[] generatePath()  {
//         gpos.resetStartProcess();
//         gneg.resetStartProcess();
//         double s;
//         for (int i=1; i<=d; i++) {
//            s = x0 + gpos.nextObservation() - gneg.nextObservation();
//            path[gpos.getCurrentObservationIndex()] = s;
//            // Note: we must get the observCounter from gpos in the case that
//            // the process is generated by a Gamma bridge
//            // The observCounter of gneg should be the same because both
//            // gamma processes should be generated the same way
//         }
//         observationIndex   = d;
//         observationCounter = d;
//         return path;
//     }

/**
 * Generates, returns and saves the path. To do so, the path of
 * @f$\Gamma^+@f$ is first generated and then the path of @f$\Gamma^-@f$.
 * This is not the optimal way of proceeding in order to reduce the variance
 * in QMC simulations; for that, use `generatePath(double[] uniform01)`
 * instead.
 */
public double[] generatePath() {
        double[] pathUP = gpos.generatePath();
        double[] pathDOWN = gneg.generatePath();

        for (int i=0; i<d; i++) {
           path[i+1] = x0 + pathUP[i+1] - pathDOWN[i+1];
        }
        observationIndex   = d;
        observationCounter = d;
        return path;
    }

   /**
    * Similar to the usual `generatePath()`, but here the uniform random
    * numbers used for the simulation must be provided to the method. This
    * allows to properly use the uniform random variates in QMC
    * simulations. This method divides the table of uniform random numbers
    * `uniform01` in two smaller tables, the first one containing the odd
    * indices of `uniform01` are used to generate the path of
    * @f$\Gamma^+@f$ and the even indices are used to generate the path
    * of @f$\Gamma^-@f$. This way of proceeding further reduces the
    * variance for QMC simulations.
    */
   public double[] generatePath (double[] uniform01) {
        int dd = uniform01.length;
        int d = dd / 2;

        if (dd % 2 != 0) {
           throw new IllegalArgumentException (
                     "The Array uniform01 must have a even length");
        }

        double[] QMCpointsUP = new double[d];
        double[] QMCpointsDW = new double[d];

        for(int i = 0; i < d; i++){
            QMCpointsUP[i] = uniform01[2*i];  // keeps the odd numbers for the gamma process
            QMCpointsDW[i] = uniform01[2*i + 1]; // and the even for the BM process
        }
        gpos.resetStartProcess();
        gneg.resetStartProcess();

        double[] pathUP = gpos.generatePath(QMCpointsUP);
        double[] pathDOWN = gneg.generatePath(QMCpointsDW);

        for (int i=0; i<d; i++) {
           path[i+1] = x0 + pathUP[i+1] - pathDOWN[i+1];
        }
        observationIndex   = d;
        observationCounter = d;
        return path;
    }

   /**
    * Sets the observation times on the `VarianceGammaProcessDiff` as
    * usual, but also applies the `resetStartProcess` method to the two
    * @ref GammaProcess objects used to generate this process.
    */
   public void resetStartProcess() {
        observationIndex   = 0;
        observationCounter = 0;
        gpos.resetStartProcess();
        gneg.resetStartProcess();
    }

   /**
    * Returns a reference to the  @ref GammaProcess object `gpos` used to
    * generate the @f$\Gamma^+@f$ component of the process.
    */
   public GammaProcess getGpos() {
        return gpos;
    }

   /**
    * Returns a reference to the  @ref GammaProcess object `gneg` used to
    * generate the @f$\Gamma^-@f$ component of the process.
    */
   public GammaProcess getGneg() {
        return gneg;
    }


    protected void init() {
        // super.init() is not called because the init() in VarianceGammaProcess
        // is very different.
        mup = 0.5 * (Math.sqrt (theta*theta + 2*sigma*sigma/nu) + theta);
        mun = 0.5 * (Math.sqrt (theta*theta + 2*sigma*sigma/nu) - theta);
        nup = mup * mup * nu;
        nun = mun * mun * nu;
        if (observationTimesSet) {
            path[0] = x0;
            gpos.setParams(t[0], mup, nup);
            gneg.setParams(t[0], mun, nun);
        }
    }

/**
 * Sets the observation times on the `VarianceGammaProcesDiff` as usual, but
 * also sets the observation times of the underlying  @ref GammaProcess ’es.
 */
public void setObservationTimes (double t[], int d) {
         gpos.setObservationTimes(t, d);
         gneg.setObservationTimes(t, d);
         super.setObservationTimes(t, d);
// the initial value is set to t[0] in the init, which is called in
// super.setObservationTimes(t, d).
     }

   /**
    * Returns the `RandomStream` of the @f$\Gamma^+@f$ process.
    */
   public RandomStream getStream() {
      return gpos.getStream ();
   }

   /**
    * Sets the  @ref umontreal.ssj.rng.RandomStream of the two
    * @ref GammaProcess ’es to `stream`.
    */
   public void setStream (RandomStream stream) {
         gpos.setStream(stream);
         gneg.setStream(stream);
   }

}