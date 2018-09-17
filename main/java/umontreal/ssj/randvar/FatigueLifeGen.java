/*
 * Class:        FatigueLifeGen
 * Description:  random variate generators for the fatigue life distribution 
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
 * This class implements random variate generators for the *fatigue life*
 * distribution @cite tBIR69a&thinsp; with location parameter @f$\mu@f$,
 * scale parameter @f$\beta@f$ and shape parameter @f$\gamma@f$. The
 * density function of this distribution is
 * @anchor REF_randvar_FatigueLifeGen_eq_fFatigueLife
 * @f[
 *   f(x) = \left[\frac{\sqrt{(x - \mu)/{\beta}} + \sqrt{{\beta}/{(x - \mu)}}}{2\gamma(x - \mu)}\right] \phi\left(\frac{\sqrt{{(x - \mu)}/{\beta}} - \sqrt{{\beta}/{(x - \mu)}}}{\gamma}\right), \qquad x > \mu\tag{fFatigueLife}
 * @f]
 * where @f$\phi@f$ is the probability density of the standard normal
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class FatigueLifeGen extends RandomVariateGen {
   protected double mu;
   protected double beta;
   protected double gamma;

   /**
    * Creates a *fatigue life* random variate generator with parameters
    * @f$\mu=@f$ `mu`, @f$\beta@f$ = `beta` and @f$\gamma@f$ = `gamma`,
    * using stream `s`.
    */
   public FatigueLifeGen (RandomStream s, double mu, double beta,
                                          double gamma) {
      super (s, new FatigueLifeDist(mu, beta, gamma));
      setParams (mu, beta, gamma);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public FatigueLifeGen (RandomStream s, FatigueLifeDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getBeta(), dist.getGamma());
   }

   /**
    * Generates a variate from the *fatigue life* distribution with
    * location parameter @f$\mu@f$, scale parameter @f$\beta@f$ and
    * shape parameter @f$\gamma@f$.
    */
   public static double nextDouble (RandomStream s, double mu, double beta,
                                    double gamma) {
      return FatigueLifeDist.inverseF (mu, beta, gamma, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Returns the parameter @f$\gamma@f$ of this object.
    */
   public double getGamma() {
      return gamma;
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Sets the parameters @f$\mu@f$, @f$\beta@f$ and @f$\gamma@f$ of
    * this object.
    */
   protected void setParams (double mu, double beta, double gamma) {
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (gamma <= 0.0)
         throw new IllegalArgumentException ("gamma <= 0");
      
      this.mu = mu;
      this.beta = beta;
      this.gamma = gamma;
   }
}