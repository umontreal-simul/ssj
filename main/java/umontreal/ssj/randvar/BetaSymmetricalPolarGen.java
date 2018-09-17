/*
 * Class:        BetaSymmetricalPolarGen
 * Description:  symmetrical beta random variate generators using
                 Ulrich's polar method 
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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
import umontreal.ssj.probdist.BetaSymmetricalDist;

/**
 * This class implements *symmetrical beta* random variate generators using
 * Ulrich’s polar method @cite rULR84a&thinsp;. The method generates two
 * uniform random variables @f$x \in[0, 1]@f$ and @f$y \in[-1, 1]@f$ until
 * @f$x^2 + y^2 \le1@f$. Then it returns
 * @anchor REF_randvar_BetaSymmetricalPolarGen_eq_beta_ulrich
 * @f[
 *   \frac{1}{2} + \frac{xy}{S}\sqrt{1 - S^{2/(2\alpha- 1)}} \tag{eq.beta.ulrich}
 * @f]
 * where @f$S = x^2 + y^2@f$, and @f$\alpha@f$ is the shape parameter of the
 * beta distribution. The method is valid only when @f$\alpha> 1/2@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class BetaSymmetricalPolarGen extends BetaSymmetricalGen {
   private double afactor;      // = 2/(2*alpha - 1)
   private RandomStream stream2;

   /**
    * Creates a symmetrical beta random variate generator with parameter
    * @f$\alpha=@f$ `alpha`, using stream `s1` to generate @f$x@f$ and
    * stream `s2` to generate @f$y@f$, as in (
    * {@link REF_randvar_BetaSymmetricalPolarGen_eq_beta_ulrich
    * eq.beta.ulrich} ) above. Restriction: @f$\alpha> 1/2@f$.
    */
   public BetaSymmetricalPolarGen (RandomStream s1, RandomStream s2,
                                   double alpha) {
      super (s1, null);
      stream2 = s2;
      if (alpha <= 0.5)
         throw new IllegalArgumentException ("  must have alpha > 1/2");
      afactor = 2.0/(2.0*alpha - 1.0);
      setParams (alpha, alpha, 0.0, 1.0);
   }

   /**
    * Creates a symmetrical beta random variate generator with parameter
    * @f$\alpha=@f$ `alpha`, using stream `s1` to generate @f$x@f$ and
    * @f$y@f$, as in (
    * {@link REF_randvar_BetaSymmetricalPolarGen_eq_beta_ulrich
    * eq.beta.ulrich} ) above. Restriction: @f$\alpha> 1/2@f$.
    */
   public BetaSymmetricalPolarGen (RandomStream s1, double alpha) {
      this (s1, s1, alpha);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s1` to generate @f$x@f$ and stream `s2` to generate @f$y@f$, as in
    * (
    * {@link REF_randvar_BetaSymmetricalPolarGen_eq_beta_ulrich
    * eq.beta.ulrich} ) above. Restriction: `dist` must have @f$\alpha>
    * 1/2@f$.
    */
   public BetaSymmetricalPolarGen (RandomStream s1, RandomStream s2,
                                   BetaSymmetricalDist dist) {
      super (s1, dist);
      stream2 = s2;
      double alp = dist.getAlpha();
      if (alp <= 0.5)
         throw new IllegalArgumentException ("  must have alpha > 1/2");
      afactor = 2.0/(2.0*dist.getAlpha() - 1.0);
      setParams (alp, alp, 0.0, 1.0);
   }

   /**
    * Creates a new generator for the distribution `dist`, using only one
    * stream `s1`. Restriction: `dist` must have @f$\alpha> 1/2@f$.
    */
   public BetaSymmetricalPolarGen (RandomStream s1,
                                   BetaSymmetricalDist dist) {
      this (s1, s1, dist);
   }

   /**
    * Generates a random number using Ulrich’s polar method. Stream `s1`
    * generates @f$x@f$ and stream `s2` generates @f$y@f$ [see eq. (
    * {@link REF_randvar_BetaSymmetricalPolarGen_eq_beta_ulrich
    * eq.beta.ulrich} )]. Restriction: @f$\alpha> 1/2@f$.
    */
   public static double nextDouble (RandomStream s1, RandomStream s2,
                                    double alpha)  {
      double u, v, S;
      do {
         u = s1.nextDouble();
         v = -1.0 + 2.0*s2.nextDouble();
         S = u*u + v*v;
      } while (S > 1.0);
      return 0.5 + u*v/S* Math.sqrt(1.0 - Math.pow(S, 2.0/(2.0*alpha - 1.0)));
   }

   /**
    * Generates a random number by Ulrich’s polar method using stream `s`.
    * Restriction: @f$\alpha> 1/2@f$.
    */
   public static double nextDouble (RandomStream s, double alpha)  {
      return nextDouble (s, s, alpha);
   }


   public double nextDouble() {
      // Generates a random number using Ulrich's polar method.
      double u, v, S;
      do {
         u = stream.nextDouble();
         v = -1.0 + 2.0*stream2.nextDouble();
         S = u*u + v*v;
      } while (S > 1.0);
      return 0.5 + u*v/S* Math.sqrt(1.0 - Math.pow(S, afactor));
  }

/**
 * Returns stream `s2` associated with this object.
 */
public RandomStream getStream2() {
      return stream2;
   }

}