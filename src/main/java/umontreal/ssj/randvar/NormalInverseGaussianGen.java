/*
 * Class:        NormalInverseGaussianGen
 * Description:  random variate generators for the normal inverse gaussian distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        June 2008
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
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements random variate generators for the *normal inverse
 * gaussian* (@f$\mathcal{NIG}@f$) distribution. See the definition of
 * @ref umontreal.ssj.probdist.NormalInverseGaussianDist
 * in package `probdist`.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalInverseGaussianGen extends RandomVariateGen {
   protected double mu;
   protected double delta = -1.0;
   protected double alpha = -1.0;
   protected double beta = -2.0;
   protected double gamma = -1.0;

   /**
    * Creates an *normal inverse gaussian* random variate generator with
    * parameters @f$\alpha@f$ = `alpha`, @f$\beta@f$ = `beta`,
    * @f$\mu@f$ = `mu` and @f$\delta@f$ = `delta`, using stream `s`.
    */
   public NormalInverseGaussianGen (RandomStream s, double alpha,
                                    double beta, double mu, double delta) {
      super (s, new NormalInverseGaussianDist(alpha, beta, mu, delta));
      setParams (alpha, beta, mu, delta);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public NormalInverseGaussianGen (RandomStream s,
                                    NormalInverseGaussianDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getAlpha(), dist.getBeta(), dist.getMu(),
                    dist.getDelta());
   }

   /**
    * NOT IMPLEMENTED. Use the daughter classes.
    */
   public static double nextDouble (RandomStream s, double alpha,
                                    double beta, double mu, double delta) {
      return NormalInverseGaussianDist.inverseF (alpha, beta, mu, delta,
                                                 s.nextDouble());
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\delta@f$ of this object.
    */
   public double getDelta() {
      return delta;
   }

   /**
    * Sets the parameters @f$\alpha@f$, @f$\beta@f$, @f$\mu@f$ and
    * @f$\delta@f$ of this object.
    */
   public void setParams (double alpha, double beta, double mu,
                          double delta) {
      if (delta <= 0.0)
         throw new IllegalArgumentException ("delta <= 0");
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (Math.abs(beta) >= alpha)
         throw new IllegalArgumentException ("|beta| >= alpha");

      gamma = Math.sqrt(alpha*alpha - beta*beta);

      this.mu = mu;
      this.delta = delta;
      this.beta = beta;
      this.alpha = alpha;
   }
 
}