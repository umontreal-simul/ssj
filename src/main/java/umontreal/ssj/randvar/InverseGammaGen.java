/*
 * Class:        InverseGammaGen
 * Description:  random variate generators for the inverse gamma distribution 
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
 * This class implements random variate generators for the *inverse gamma*
 * distribution with shape parameter @f$\alpha> 0@f$ and scale parameter
 * @f$\beta> 0@f$. The density function of this distribution is
 * @anchor REF_randvar_InverseGammaGen_eq_finvgam
 * @f[
 *   f(x) = \left\{\begin{array}{ll}
 *    \displaystyle\frac{\beta^{\alpha}}{ \Gamma(\alpha)} \frac{ e^{-\beta/ x}}{ x^{\alpha+ 1}} 
 *    & 
 *    \quad\mbox{for } x > 0 
 *    \\ 
 *    0 
 *    & 
 *    \quad\mbox{otherwise,} 
 *   \end{array} \right. \tag{finvgam}
 * @f]
 * where @f$\Gamma@f$ is the gamma function.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class InverseGammaGen extends RandomVariateGen {
   protected double alpha;
   protected double beta;

   /**
    * Creates an inverse gamma random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta`, using stream `s`.
    */
   public InverseGammaGen (RandomStream s, double alpha, double beta) {
      super (s, new InverseGammaDist(alpha, beta));
      setParams(alpha, beta);
   }

   /**
    * Creates an inverse gamma random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta= 1@f$, using stream `s`.
    */
   public InverseGammaGen (RandomStream s, double alpha) {
      this (s, alpha, 1.0);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public InverseGammaGen (RandomStream s, InverseGammaDist dist) {
      super (s, dist);
      if (dist != null)
         setParams(dist.getAlpha(), dist.getBeta());
   }

   /**
    * Generates a variate from the inverse gamma distribution with shape
    * parameter @f$\alpha> 0@f$ and scale parameter @f$\beta> 0@f$.
    */
   public static double nextDouble (RandomStream s,
                                    double alpha, double beta) {
      return InverseGammaDist.inverseF (alpha, beta, s.nextDouble());
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


   protected void setParams (double alpha, double beta) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      this.alpha = alpha;
      this.beta = beta;
   }
}