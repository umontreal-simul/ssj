/*
 * Class:        CIRProcess
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
 * This class represents a *CIR* (Cox, Ingersoll, Ross) process
 * @cite fCOX85a&thinsp; @f$\{X(t) : t \geq0 \}@f$, sampled at times @f$0 =
 * t_0 < t_1 < \cdots< t_d@f$. This process obeys the stochastic
 * differential equation
 * @anchor REF_stochprocess_CIRProcess_eq_cir
 * @f[
 *   dX(t) = \alpha(b - X(t)) dt + \sigma\sqrt{X(t)}  dB(t) \tag{cir}
 * @f]
 * with initial condition @f$X(0)= x_0@f$, where @f$\alpha@f$, @f$b@f$ and
 * @f$\sigma@f$ are positive constants, and @f$\{B(t),  t\ge0\}@f$ is a
 * standard Brownian motion (with drift 0 and volatility 1). This process is
 * *mean-reverting* in the sense that it always tends to drift toward its
 * general mean @f$b@f$. The process is generated using the sequential
 * technique @cite fGLA04a&thinsp; (p. 122)
 * @anchor REF_stochprocess_CIRProcess_eq_cir_seq
 * @f[
 *   X(t_j) = \frac{\sigma^2\left(1 - e^{-\alpha(t_j - t_{j-1})}\right)}{4\alpha} \chi^{\prime 2}_{\nu}\left(\frac{4\alpha e^{-\alpha(t_j - t_{j-1}) } X(t_{j-1})}{\sigma^2\left(1 - e^{-\alpha(t_j - t_{j-1})}\right)}\right), \tag{cir-seq}
 * @f]
 * where @f$\nu= 4b\alpha/\sigma^2@f$, and @f$\chi^{\prime
 * 2}_{\nu}(\lambda)@f$ is a noncentral chi-square random variable with
 * @f$\nu@f$ degrees of freedom and noncentrality parameter @f$\lambda@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class CIRProcess extends StochasticProcess {
    private RandomStream stream;
    protected ChiSquareNoncentralGen gen;
    protected double alpha,
                     beta,
                     sigma,
                     nu;     // number of degrees of freedom
    // Precomputed values
    protected double[] parc,
                       parlam;

   /**
    * Constructs a new `CIRProcess` with parameters @f$\alpha=
    * \mathtt{alpha}@f$, @f$b@f$, @f$\sigma= \mathtt{sigma}@f$ and
    * initial value @f$X(t_0) = \mathtt{x0}@f$. The noncentral chi-square
    * variates @f$\chi^{\prime2}_{\nu}(\lambda)@f$ will be generated
    * by inversion using the stream `stream`.
    */
   public CIRProcess (double x0, double alpha, double b, double sigma,
                      RandomStream stream) {
      this (x0, alpha, b, sigma, new ChiSquareNoncentralGen (stream, null));
    }

   /**
    * The noncentral chi-square variate generator `gen` is specified
    * directly instead of specifying the stream. `gen` can use a method
    * other than inversion.
    */
   public CIRProcess (double x0, double alpha, double b, double sigma,
                      ChiSquareNoncentralGen gen) {
      this.alpha = alpha;
      this.beta  = b;
      this.sigma = sigma;
      this.x0    = x0;
      nu = 4.0*b*alpha/(sigma*sigma);
      this.gen   = gen;
      stream = gen.getStream();
    }


   public double nextObservation() {
      double xOld = path[observationIndex];
      double lambda = xOld * parlam[observationIndex];
      double x;
      if (gen.getClass() == ChiSquareNoncentralPoisGen.class)
         x = parc[observationIndex] *
             ChiSquareNoncentralPoisGen.nextDouble(stream, nu, lambda);
      else if (gen.getClass() == ChiSquareNoncentralGamGen.class)
         x = parc[observationIndex] *
             ChiSquareNoncentralGamGen.nextDouble(stream, nu, lambda);
      else
         x = parc[observationIndex] *
             ChiSquareNoncentralGen.nextDouble(stream, nu, lambda);
      observationIndex++;
      path[observationIndex] = x;
      return x;
   }

/**
 * Generates and returns the next observation at time @f$t_{j+1} =
 * \mathtt{nextTime}@f$, using the previous observation time @f$t_j@f$
 * defined earlier (either by this method or by
 * <tt>setObservationTimes</tt>), as well as the value of the previous
 * observation @f$X(t_j)@f$. *Warning*: This method will reset the
 * observations time @f$t_{j+1}@f$ for this process to `nextTime`. The user
 * must make sure that the @f$t_{j+1}@f$ supplied is @f$\geq t_j@f$.
 */
public double nextObservation (double nextTime) {
      double previousTime = t[observationIndex];
      double xOld = path[observationIndex];
      observationIndex++;
      t[observationIndex] = nextTime;
      double dt = nextTime - previousTime;
      double x = nextObservation (xOld, dt);
      path[observationIndex] = x;
      return x;
    }

   /**
    * Generates an observation of the process in `dt` time units, assuming
    * that the process has value @f$x@f$ at the current time. Uses the
    * process parameters specified in the constructor. Note that this
    * method does not affect the sample path of the process stored
    * internally (if any).
    */
   public double nextObservation (double x, double dt) {
      double c = -Math.expm1(-alpha * dt) * sigma * sigma / (4.0*alpha);
      double lambda = x * Math.exp(-alpha * dt) / c;
      if (gen.getClass() == ChiSquareNoncentralPoisGen.class)
         x = c * ChiSquareNoncentralPoisGen.nextDouble(stream, nu, lambda);
      else if (gen.getClass() == ChiSquareNoncentralGamGen.class)
         x = c * ChiSquareNoncentralGamGen.nextDouble(stream, nu, lambda);
      else
         x = c * ChiSquareNoncentralGen.nextDouble(stream, nu, lambda);
      return x;
    }
public double[] generatePath() {
      double xOld = x0;
      double x, lambda;
      int j;

      if (gen.getClass() == ChiSquareNoncentralPoisGen.class) {
         for (j = 0; j < d; j++) {
            lambda = xOld * parlam[j];
            x = parc[j] * ChiSquareNoncentralPoisGen.nextDouble(stream, nu, lambda);
            path[j + 1] = x;
            xOld = x;
         }

      } else if (gen.getClass() == ChiSquareNoncentralGamGen.class) {
         for (j = 0; j < d; j++) {
            lambda = xOld * parlam[j];
            x = parc[j] * ChiSquareNoncentralGamGen.nextDouble(stream, nu, lambda);
            path[j + 1] = x;
            xOld = x;
         }

      } else {
         for (j = 0; j < d; j++) {
            lambda = xOld * parlam[j];
            x = parc[j] * ChiSquareNoncentralGen.nextDouble(stream, nu, lambda);
            path[j + 1] = x;
            xOld = x;
         }
      }

      observationIndex = d;
      return path;
   }

/**
 * Generates a sample path of the process at all observation times, which are
 * provided in array `t`. Note that `t[0]` should be the observation time of
 * `x0`, the initial value of the process, and `t[]` should have at least
 * @f$d+1@f$ elements (see the `setObservationTimes` method).
 */
public double[] generatePath (RandomStream stream) {
      gen.setStream (stream);
      return generatePath();
   }

   /**
    * Resets the parameters @f$X(t_0) = \mathtt{x0}@f$, @f$\alpha=
    * \mathtt{alpha}@f$, @f$b = \mathtt{b}@f$ and @f$\sigma=
    * \mathtt{sigma}@f$ of the process. *Warning*: This method will
    * recompute some quantities stored internally, which may be slow if
    * called too frequently.
    */
   public void setParams (double x0, double alpha, double b, double sigma) {
      this.alpha = alpha;
      this.beta  = b;
      this.sigma = sigma;
      this.x0    = x0;
      nu = 4.0*b*alpha/(sigma*sigma);
      if (observationTimesSet) init(); // Otherwise not needed.
    }

   /**
    * Resets the random stream of the noncentral chi-square generator to
    * `stream`.
    */
   public void setStream (RandomStream stream) { gen.setStream (stream); }

   /**
    * Returns the random stream of the noncentral chi-square generator.
    */
   public RandomStream getStream() { return gen.getStream (); }

   /**
    * Returns the value of @f$\alpha@f$.
    */
   public double getAlpha() { return alpha; }

   /**
    * Returns the value of @f$b@f$.
    */
   public double getB() { return beta; }

   /**
    * Returns the value of @f$\sigma@f$.
    */
   public double getSigma() { return sigma; }

   /**
    * Returns the noncentral chi-square random variate generator used. The
    * `RandomStream` used for that generator can be changed via
    * `getGen().setStream(stream)`, for example.
    */
   public ChiSquareNoncentralGen getGen() { return gen; }


   protected void initArrays(int d) {
      double dt, c;
      for (int j = 0; j < d; j++) {
         dt = t[j+1] - t[j];
         c = -Math.expm1(-alpha * dt)*sigma*sigma/(4.0*alpha);
         parc[j] = c;
         parlam[j] = Math.exp(-alpha * dt) / c;
      }
   }

   // This is called by setObservationTimes to precompute constants
   // in order to speed up the path generation.
   protected void init() {
       super.init();
       parc = new double[d];
       parlam = new double[d];
       initArrays(d);
    }

}