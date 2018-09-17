/*
 * Class:        ParetoGen
 * Description:  random variate generators for the Pareto distribution
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
 * This class implements random variate generators for one of the *Pareto*
 * distributions, with parameters @f$\alpha>0@f$ and @f$\beta>0@f$. Its
 * density function is
 * @anchor REF_randvar_ParetoGen_eq_fpareto
 * @f[
 *   f (x) = \left\{\begin{array}{ll}
 *    {\displaystyle\frac{\alpha\beta^{\alpha}}{x^{\alpha+1}}} 
 *    & 
 *   \mbox{ for }x>\beta
 *    \\ 
 *    0 
 *    & 
 *    \mbox{ for }x\le\beta
 *   \end{array}\right. \tag{fpareto}
 * @f]
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ParetoGen extends RandomVariateGen {
   protected double alpha;
   protected double beta;

   /**
    * Creates a Pareto random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta= @f$ `beta`, using stream `s`.
    */
   public ParetoGen (RandomStream s, double alpha, double beta) {
      super (s, new ParetoDist(alpha, beta));
      setParams(alpha, beta);
   }

   /**
    * Creates a Pareto random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta= 1@f$, using stream `s`.
    */
   public ParetoGen (RandomStream s, double alpha) {
      this (s, alpha, 1.0);
   }

   /**
    * Creates a new generator for the Pareto distribution `dist` and
    * stream `s`.
    */
   public ParetoGen (RandomStream s, ParetoDist dist) {
      super (s, dist);
      if (dist != null)
         setParams(dist.getAlpha(), dist.getBeta());
   }

   /**
    * Generates a new variate from the Pareto distribution with parameters
    * @f$\alpha= @f$&nbsp;`alpha` and @f$\beta= @f$&nbsp;`beta`, using
    * stream `s`.
    */
   public static double nextDouble (RandomStream s,
                                    double alpha, double beta) {
      return ParetoDist.inverseF (alpha, beta, s.nextDouble());
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