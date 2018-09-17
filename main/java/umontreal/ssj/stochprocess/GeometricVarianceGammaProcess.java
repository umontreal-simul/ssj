/*
 * Class:        GeometricVarianceGammaProcess
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Pierre Tremblay
 * @since        July 2003
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
 * This class represents a *geometric variance gamma* process @f$S(t)@f$ (see
 * @cite fMAD98a&thinsp; (page 86)). This stochastic process is defined by
 * the equation
 * @anchor REF_stochprocess_GeometricVarianceGammaProcess_GeoVGeqn
 * @f[
 *   S(t) = S(0) \mbox{ exp}(\mu t + X(t; \sigma, \nu, \theta) + \omega t), \tag{GeoVGeqn}
 * @f]
 * where @f$X@f$ is a variance gamma process and
 * @anchor REF_stochprocess_GeometricVarianceGammaProcess_omegaEqn
 * @f[
 *   \omega= (1/\nu) \mbox{ ln}( 1 - \theta\nu- \sigma^2 \nu/2). \tag{omegaEqn}
 * @f]
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class GeometricVarianceGammaProcess extends StochasticProcess {
    protected VarianceGammaProcess vargamma;
    protected double        theta,
                            nu,
                            mu,
                            sigma,
                            omega,
                            muPlusOmega;
    protected double[]      mudt;

   /**
    * Constructs a new `GeometricVarianceGammaProcess` with parameters
    * @f$\theta= \mathtt{theta}@f$, @f$\sigma= \mathtt{sigma}@f$,
    * @f$\nu= \mathtt{nu}@f$, @f$\mu= \mathtt{mu}@f$ and initial value
    * @f$S(t_0) = \mathtt{s0}@f$. The `stream` is used to generate the
    * @ref VarianceGammaProcess object used to implement @f$X@f$ in (
    * {@link REF_stochprocess_GeometricVarianceGammaProcess_GeoVGeqn
    * GeoVGeqn} ).
    */
   public GeometricVarianceGammaProcess (double s0, double theta,
                                         double sigma, double nu,
                                         double mu, RandomStream stream) {
        vargamma = new VarianceGammaProcess (0.0, theta, sigma, nu, stream);
        setParams (s0, theta, sigma, nu, mu);
    }

   /**
    * Constructs a new `GeometricVarianceGammaProcess`. The parameters
    * @f$\theta, \sigma, \nu@f$ are set to the parameters of the
    * @ref VarianceGammaProcess `vargamma`. The parameter @f$\mu@f$ is
    * set to `mu` and the initial values @f$S(t_0) = \mathtt{s0}@f$.
    */
   public GeometricVarianceGammaProcess (double s0, double mu,
                                         VarianceGammaProcess vargamma) {
        this.vargamma = vargamma;
        setParams (s0, vargamma.getTheta (), vargamma.getSigma (),
                   vargamma.getNu (), mu);
   }
   
   public double nextObservation() {
        double nextX  = vargamma.nextObservation();
        observationIndex = vargamma.getCurrentObservationIndex();
        // Could be different than simply 'observationIndex++' because of the
        // possibility of Gamma/Brownian bridge
        observationCounter++;

        double s = x0 * Math.exp (muPlusOmega * (t[observationIndex] - t[0])
                                  + nextX);
        path[observationIndex] = s;
        return s;
    }

   public double[] generatePath() {
        double s = x0;
        resetStartProcess();
        double[] vgpath = vargamma.generatePath();
        for (int i = 0; i < d; i++) {
            s *= Math.exp (mudt[i] + vgpath[i+1] - vgpath[i]);
            path[i+1] = s;
        }
        observationIndex = d;
        observationCounter++;
        return path;
    }
    // allows the user to create a path by specifiying the uniform random numbers to be used
   public double[] generatePath (double[] uniform01) {
        double s = x0;
        resetStartProcess();

        double[] vgpath = vargamma.generatePath(uniform01);
        for (int i = 0; i < d; i++) {
            s *= Math.exp (mudt[i] + vgpath[i+1] - vgpath[i]);
            path[i+1] = s;
        }
        observationIndex = d;
        observationCounter++;
        return path;
    }


    // method not verified by JS...  old stuff
   public double getCurrentUpperBound()  {
        // Find index for last observation generated (chronologically)
        int j = 0; // By default, t0 !
        int i = observationCounter - 1;
        double tForIthObserv;
        while (i > 0) {
            tForIthObserv = t[observationIndexFromCounter[i]];
            if (tForIthObserv <= t[observationCounter] && tForIthObserv > t[j])
                j = i;
            i--;
        }

        // Calculate bound following recipe
        double u = 0.0;
        GammaProcess gpos = ((VarianceGammaProcessDiff) vargamma).getGpos();
        double[] gposPath = gpos.getPath();
        double deltaGpos = gposPath[observationIndex] - gposPath[j];
        double s = path[observationIndex];
        if (muPlusOmega < 0)
             u = s * Math.exp (deltaGpos);
        else u = s * Math.exp (muPlusOmega * (t[observationIndex] - t[j])
                               + deltaGpos);
        return u;
    }

   /**
    * Resets the `GeometricaVarianceGammaProcess`, but also applies the
    * `resetStartProcess` method to the  @ref VarianceGammaProcess object
    * used to generate this process.
    */
   public void resetStartProcess() {
        observationIndex   = 0;
        observationCounter = 0;
        vargamma.resetStartProcess();
    }

   /**
    * Sets the parameters @f$S(t_0) = \mathtt{s0}@f$, @f$\theta=
    * \mathtt{theta}@f$, @f$\sigma= \mathtt{sigma}@f$, @f$\nu=
    * \mathtt{nu}@f$ and @f$\mu= \mathtt{mu}@f$ of the process.
    * *Warning*: This method will recompute some quantities stored
    * internally, which may be slow if called repeatedly.
    */
   public void setParams (double s0, double theta, double sigma, double nu,
                          double mu) {
        this.x0    = s0;
        this.theta = theta;
        this.sigma = sigma;
        this.nu    = nu;
        this.mu    = mu;
        if (observationTimesSet) init(); // Otherwise no need to.
    }

   /**
    * Returns the value of the parameter @f$\theta@f$.
    */
   public double getTheta() { return theta; }

   /**
    * Returns the value of the parameter @f$\mu@f$.
    */
   public double getMu() { return mu; }

   /**
    * Returns the value of the parameter @f$\nu@f$.
    */
   public double getNu() { return nu; }

   /**
    * Returns the value of the parameter @f$\sigma@f$.
    */
   public double getSigma() { return sigma; }

   /**
    * Returns the value of the quantity @f$\omega@f$ defined in (
    * {@link REF_stochprocess_GeometricVarianceGammaProcess_omegaEqn
    * omegaEqn} ).
    */
   public double getOmega() { return omega; }

   /**
    * Returns a reference to the variance gamma process @f$X@f$ defined in
    * the constructor.
    */
   public VarianceGammaProcess getVarianceGammaProcess() {
        return vargamma;
}


    protected void init() {
        super.init();
        if (1 <= theta*nu + sigma*sigma*nu / 2.0)
           throw new IllegalArgumentException ("theta*nu + sigma*sigma*nu / 2 >= 1");
        omega = Math.log (1 - theta*nu - sigma*sigma*nu / 2.0) / nu;
        muPlusOmega = mu + omega;

        if (observationTimesSet) {
            // Telling the variance gamma proc. about the observ. times
            vargamma.setObservationTimes (t, d);

            // We need to know in which order the observations are generated
            this.observationIndexFromCounter
                = vargamma.getArrayMappingCounterToIndex();

            mudt = new double[d];
            for (int i = 0; i < d; i++)
                mudt[i] = muPlusOmega * (t[i+1] - t[i]);
        }
    }


   public void setStream (RandomStream stream)  {
        vargamma.setStream(stream);
   }


   public RandomStream getStream () {
      return vargamma.getStream();
   }

}