/*
 * Class:        GeometricBrownianMotion
 * Description:
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
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;

/**
 * Represents a *geometric Brownian motion* (GBM) process @f$\{S(t),
 * t\ge0\}@f$, which evolves according to the stochastic differential
 * equation
 * @anchor REF_stochprocess_GeometricBrownianMotion_eq_GBM
 * @f[
 *   dS(t) = \mu S(t) dt + \sigma S(t) dB(t), \tag{GBM}
 * @f]
 * where @f$\mu@f$ and @f$\sigma@f$ are the drift and volatility
 * parameters, and @f$\{B(t),  t\ge0\}@f$ is a standard Brownian motion (for
 * which @f$B(t)\sim N(0,t)@f$). This process can also be written as the
 * exponential of a Brownian motion:
 * @anchor REF_stochprocess_GeometricBrownianMotion_eq_GBM2
 * @f[
 *   S(t) = S(0) \exp\left[ (\mu- \sigma^2/2) t + \sigma B(t) \right] = S(0) \exp\left[ X(t) \right], \tag{GBM2}
 * @f]
 * where @f$X(t) = (\mu- \sigma^2/2) t + \sigma B(t)@f$. The GBM process
 * is simulated by simulating the BM process @f$X@f$ and taking the
 * exponential. This BM process is stored internally.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class GeometricBrownianMotion extends StochasticProcess {

    protected NormalGen      gen;
    protected BrownianMotion bm;   // The underlying BM process X.
    protected double         mu,
                             sigma;
    protected double[]       mudt;

   /**
    * Same as `GeometricBrownianMotion (s0, mu, sigma, new BrownianMotion
    * (0.0, 0.0, 1.0, stream))`.
    */
   public GeometricBrownianMotion (double s0, double mu, double sigma,
                                   RandomStream stream) {
        this (s0, mu, sigma, new BrownianMotion (0.0, 0.0, 1.0, stream));
    }

   /**
    * Constructs a new `GeometricBrownianMotion` with parameters @f$\mu=
    * \mathtt{mu}@f$, @f$\sigma= \mathtt{sigma}@f$, and @f$S(t_0) =
    * \mathtt{s0}@f$, using `bm` as the underlying  @ref BrownianMotion.
    * The parameters of `bm` are automatically reset to
    * @f$\mu-\sigma^2/2@f$ and @f$\sigma@f$, regardless of the original
    * parameters of `bm`. The observation times are the same as those of
    * `bm`. The generation method depends on that of `bm` (sequential,
    * bridge sampling, PCA, etc.).
    */
   public GeometricBrownianMotion (double s0, double mu, double sigma,
                                   BrownianMotion bm) {
        this.bm = bm;
        setParams (s0, mu, sigma);
    }

   public void setObservationTimes (double[] t, int d) {
        this.d = d;
        super.setObservationTimes (t, d);
        bm.setObservationTimes (t, d);
    }

   public double nextObservation() {
        // Note : this implementation is general, to deal with
        // the possibility of generating bm with bridge sampling, for example.  ???

        double s = x0 * Math.exp (bm.nextObservation());
        observationIndex = bm.getCurrentObservationIndex();
        path[observationIndex] = s;
        // Could be different than simply 'observationCounter++' because of the
        // possibility of Brownian bridge

        return s;
    }

   public double[] generatePath() {
        path[0] = x0;
        bm.generatePath ();
        for (int i = 1; i <= d; ++i)
            path[i] = x0 * Math.exp (bm.getObservation(i));
        observationCounter = d;
        return path;
    }

   public double[] generatePath (RandomStream stream) {
        setStream (stream);
        return generatePath();
    }

/**
 * Same as in `StochasticProcess`, but also invokes `resetStartProcess` for
 * the underlying `BrownianMotion` object.
 */
public void resetStartProcess() {
        observationCounter = 0;
        bm.resetStartProcess();
    }

   /**
    * Sets the parameters @f$S(t_0) = \mathtt{s0}@f$, @f$\mu=
    * \mathtt{mu}@f$ and @f$\sigma= \mathtt{sigma}@f$ of the process.
    * *Warning*: This method will recompute some quantities stored
    * internally, which may be slow if called repeatedly.
    */
   public void setParams (double s0, double mu, double sigma) {
        this.x0    = s0;
        this.mu    = mu;
        this.sigma = sigma;
        bm.setParams (0.0, mu - 0.5 * sigma * sigma, sigma);
        if (observationTimesSet) init(); // Otherwise not needed.
    }

   /**
    * Resets the  @ref umontreal.ssj.rng.RandomStream for the underlying
    * Brownian motion to `stream`.
    */
   public void setStream (RandomStream stream) { (bm.gen).setStream (stream); }

   /**
    * Returns the  @ref umontreal.ssj.rng.RandomStream for the underlying
    * Brownian motion.
    */
   public RandomStream getStream() { return (bm.gen).getStream (); }

   /**
    * Returns the value of @f$\mu@f$.
    */
   public double getMu() { return mu; }

   /**
    * Returns the value of @f$\sigma@f$.
    */
   public double getSigma() { return sigma; }

   /**
    * Returns the  @ref umontreal.ssj.randvar.NormalGen used.
    */
   public NormalGen getGen() { return gen; }

   /**
    * Returns a reference to the  @ref BrownianMotion object used to
    * generate the process.
    */
   public BrownianMotion getBrownianMotion() {
        return bm;
    }


    protected void init() {
        super.init();   // Maybe useless...
    }

}
