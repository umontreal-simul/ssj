/*
 * Class:        HalfNormalGen
 * Description:  generator of random variates from the half-normal distribution
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
 * This class implements methods for generating random variates from the
 * *half-normal* distribution with parameters @f$\mu@f$ and @f$\sigma>
 * 0@f$. Its density is
 * @f{align*}{
 *    f(x) 
 *    & 
 *   =
 *    \frac{1}{\sigma}\sqrt{\frac{2}{\pi}}\; e^{-(x-\mu)^2/2\sigma^2}, \qquad\mbox{for } x \ge\mu. 
 *    \\ 
 *    f(x) 
 *    & 
 *   =
 *    0, \qquad\mbox{for } x < \mu.
 * @f}
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class HalfNormalGen extends RandomVariateGen {
    
   // Distribution parameters
   protected double mu;
   protected double sigma;

   /**
    * Creates a new *half-normal* generator with parameters @f$\mu=@f$
    * `mu` and @f$\sigma=@f$ `sigma`, using stream `s`.
    */
   public HalfNormalGen (RandomStream s, double mu, double sigma) {
      super (s, new HalfNormalDist (mu, sigma));
      setParams (mu, sigma);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public HalfNormalGen (RandomStream s, HalfNormalDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getSigma());
   }

   /**
    * Generates a variate from the *half-normal* distribution with
    * parameters @f$\mu= @f$&nbsp;`mu` and @f$\sigma= @f$&nbsp;`sigma`,
    * using stream `s`.
    *  @param s            the random stream
    *  @param mu           the parameter mu
    *  @param sigma        the parameter sigma
    *  @return Generates a variate from the *HalfNormal* distribution
    */
   public static double nextDouble (RandomStream s, double mu, double sigma) {
      return HalfNormalDist.inverseF (mu, sigma, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    *  @return the parameter mu
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\sigma@f$ of this object.
    *  @return the parameter mu
    */
   public double getSigma() {
      return sigma;
   }


   protected void setParams (double mu, double sigma) {
      if (sigma <= 0.0)
         throw new IllegalArgumentException ("sigma <= 0");
      this.mu = mu;
      this.sigma = sigma;
   }
}