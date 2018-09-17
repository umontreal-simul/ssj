/*
 * Class:        InverseGaussianGen
 * Description:  random variate generators for the inverse Gaussian distribution
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
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements random variate generators for the *inverse Gaussian*
 * distribution with location parameter @f$\mu> 0@f$ and scale parameter
 * @f$\lambda> 0@f$. The density function of this distribution is
 * @anchor REF_randvar_InverseGaussianGen_eq_fInverseGaussian
 * @f[
 *   f(x) = \sqrt{\frac{\lambda}{2\pi x^3}}\; e^{{-\lambda(x - \mu)^2}/{(2\mu^2x)}} \qquad\mbox{for } x > 0. \tag{fInverseGaussian}
 * @f]
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class InverseGaussianGen extends RandomVariateGen {
   protected double mu = -1.0;
   protected double lambda = -1.0;

   /**
    * Creates an *inverse Gaussian* random variate generator with
    * parameters @f$\mu= @f$ `mu` and @f$\lambda= @f$ `lambda`, using
    * stream `s`.
    */
   public InverseGaussianGen (RandomStream s, double mu, double lambda) {
      super (s, new InverseGaussianDist(mu, lambda));
      setParams (mu, lambda);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public InverseGaussianGen (RandomStream s, InverseGaussianDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getLambda());
   }

   /**
    * Generates a variate from the inverse gaussian distribution with
    * location parameter @f$\mu> 0@f$ and scale parameter @f$\lambda>
    * 0@f$.
    */
   public static double nextDouble (RandomStream s,
                                    double mu, double lambda) {
      return InverseGaussianDist.inverseF (mu, lambda, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\lambda@f$ of this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\lambda@f$ of this object.
    */
   protected void setParams (double mu, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");
      this.mu = mu;
      this.lambda = lambda;
   }
}