/*
 * Class:        BetaSymmetricalBestGen
 * Description:  symmetrical beta random variate generators using
                 Devroye's one-liner method.
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
 * Devroye’s one-liner method. It is based on Best’s relation
 * @cite rBES78a&thinsp; between a Student-@f$t@f$ variate and a symmetrical
 * beta variate:
 * @f[
 *   B_{\alpha, \alpha} \stackrel{\mathcal{L}}{=} \frac{1}{2} \left( 1 + \frac{T_{2\alpha}}{\sqrt{2\alpha+ T_{2\alpha}^2}}\right).
 * @f]
 * If @f$S@f$ is a random sign and @f$U_1@f$, @f$U_2@f$ are two independent
 * uniform @f$[0,1]@f$ random variates, then the following gives a
 * symmetrical beta variate @cite rDEV96a&thinsp;:
 * @anchor REF_randvar_BetaSymmetricalBestGen_eq_beta_best
 * @f[
 *   B_{\alpha, \alpha} \stackrel{\mathcal{L}}{=} \frac{1}{2} + \frac{S}{2 \sqrt{1 + \frac{1}{\left(U_1^{-1/\alpha} - 1\right)\cos^2(2\pi U_2)}}} \tag{eq.beta.best}
 * @f]
 * valid for any shape parameter @f$\alpha> 0@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class BetaSymmetricalBestGen extends BetaSymmetricalGen {
   private RandomStream stream2;
   private RandomStream stream3;
   private double afactor;          // = 1/alpha
   private static final double TWOPI = 2.0*Math.PI;      // = 2 Pi

   /**
    * Creates a symmetrical beta random variate generator with parameter
    * @f$\alpha=@f$ `alpha`, using stream `s1` to generate @f$U_1@f$,
    * stream `s2` to generate @f$U_2@f$ and stream `s3` to generate
    * @f$S@f$, as given in equation (
    * {@link REF_randvar_BetaSymmetricalBestGen_eq_beta_best
    * eq.beta.best} ).
    */
   public BetaSymmetricalBestGen (RandomStream s1, RandomStream s2,
                                  RandomStream s3, double alpha) {
      super (s1, null);
      stream2 = s2;
      stream3 = s3;
      afactor = 1.0/alpha;
      setParams (alpha, alpha, 0.0, 1.0);
   }

   /**
    * Creates a symmetrical beta random variate generator with parameter
    * @f$\alpha=@f$ `alpha`, using only one stream `s1` to generate
    * @f$U_1@f$, @f$U_2@f$, and @f$S@f$ as given in equation (
    * {@link REF_randvar_BetaSymmetricalBestGen_eq_beta_best
    * eq.beta.best} ).
    */
   public BetaSymmetricalBestGen (RandomStream s1, double alpha) {
     this (s1, s1, s1, alpha);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s1` to generate @f$U_1@f$, stream `s2` to generate @f$U_2@f$ and
    * stream `s3` to generate @f$S@f$ as given in equation (
    * {@link REF_randvar_BetaSymmetricalBestGen_eq_beta_best
    * eq.beta.best} ).
    */
   public BetaSymmetricalBestGen (RandomStream s1, RandomStream s2,
                                  RandomStream s3, BetaSymmetricalDist dist) {
      super (s1, dist);
      stream2 = s2;
      stream3 = s3;
      afactor = 1.0/dist.getAlpha();
      if (dist != null)
         setParams (dist.getAlpha(), dist.getAlpha(), dist.getA(), dist.getB());
   }

   /**
    * Creates a new generator for the distribution `dist`, using only one
    * stream `s1`.
    */
   public BetaSymmetricalBestGen (RandomStream s1, BetaSymmetricalDist dist) {
     this (s1, s1, s1, dist);
   }

   /**
    * Generates a random number using Devroye’s one-liner method.
    * Restriction: @f$\alpha> 0@f$.
    */
   public static double nextDouble (RandomStream s1, RandomStream s2,
                                    RandomStream s3, double alpha)  {
      double cos, temp, v, S;
      cos = Math.cos (TWOPI * s2.nextDouble());
      temp = 1.0/Math.pow(s1.nextDouble(), 1.0/alpha) - 1.0;
      v = Math.sqrt(1.0 + 1.0 / (temp*cos*cos));
      S = s3.nextDouble();
      if (S < 0.5)
         return 0.5 - 0.5/v;
      else
         return 0.5 + 0.5/v;
   }

   /**
    * Generates a random number using Devroye’s one-liner method with only
    * one stream `s`. Restriction: @f$\alpha> 0@f$.
    */
   public static double nextDouble (RandomStream s, double alpha)  {
      return nextDouble (s, s, s, alpha);
   }


   public double nextDouble() {
      // Generates a random number using Devroye's one liner method
      double cos, temp, v, S;
      cos = Math.cos (TWOPI * stream2.nextDouble());
      temp = 1.0/Math.pow(stream.nextDouble(), afactor) - 1.0;
      v = Math.sqrt(1.0 + 1.0 / (temp*cos*cos));
      S = stream3.nextDouble();
      if (S < 0.5)
         return 0.5 - 0.5/v;
      else
         return 0.5 + 0.5/v;
  }

/**
 * Returns stream `s2` associated with this object.
 */
public RandomStream getStream2() {
      return stream2;
   }

   /**
    * Returns stream `s3` associated with this object.
    */
   public RandomStream getStream3() {
      return stream3;
   }

}