/*
 * Class:        BetaGen
 * Description:  random variate generators with the beta distribution
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
 * This class implements random variate generators with the *beta*
 * distribution with shape parameters @f$\alpha> 0@f$ and @f$\beta> 0@f$,
 * over the interval @f$(a,b)@f$, where @f$a < b@f$. The density function of
 * this distribution is
 * @f[
 *   f(x) = \frac{\Gamma(\alpha+\beta)}{ \Gamma(\alpha)\Gamma(\beta)(b - a)^{\alpha+\beta-1}} (x - a)^{\alpha- 1}(b - x)^{\beta-1} \qquad\mbox{ for } a < x < b,
 * @f]
 * and @f$f(x)=0@f$ elsewhere, where @f$\Gamma(x)@f$ is the gamma function
 * defined in
 *  ( {@link REF_randvar_GammaGen_eq_Gamma Gamma} ).
 *
 * Local copies of the parameters @f$\alpha@f$, @f$\beta@f$, @f$a@f$, and
 * @f$b@f$ are maintained in this class. The (non-static) `nextDouble` method
 * simply calls `inverseF` on the distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class BetaGen extends RandomVariateGen {
    
   // Distribution parameters
   protected double p;
   protected double q;
   protected double a;
   protected double b;
   protected int gen;

   /**
    * Creates a new beta generator with parameters @f$\alpha=@f$ `alpha`
    * and @f$\beta=@f$ `beta`, over the interval
    * @f$(@f$<tt>a</tt>@f$,@f$&nbsp;<tt>b</tt>@f$)@f$, using stream `s`.
    */
   public BetaGen (RandomStream s, double alpha, double beta,
                                   double a, double b) {
      super (s, new BetaDist (alpha, beta, a, b));
      setParams (alpha, beta, a, b);
   }

   /**
    * Creates a new beta generator with parameters @f$\alpha=@f$ `alpha`
    * and @f$\beta=@f$ `beta`, over the interval @f$(0,1)@f$, using
    * stream `s`.
    */
   public BetaGen (RandomStream s, double alpha, double beta) {
      this (s, alpha, beta, 0.0, 1.0);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public BetaGen (RandomStream s, BetaDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getAlpha(), dist.getBeta(), dist.getA(), dist.getB());
   }

   /**
    * Generates a variate from the *beta* distribution with parameters
    * @f$\alpha= @f$&nbsp;`alpha`, @f$\beta= @f$&nbsp;`beta`, over the
    * interval @f$(a, b)@f$, using stream `s`.
    */
   public static double nextDouble (RandomStream s, double alpha,
                                    double beta, double a, double b) {
      return BetaDist.inverseF (alpha, beta, a, b, 15, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return p;
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return q;
   }

   /**
    * Returns the parameter @f$a@f$ of this object.
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the parameter @f$b@f$ of this object.
    */
   public double getB() {
      return b;
   }


   protected void setParams (double alpha, double beta, double aa, double bb) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (aa >= bb)
         throw new IllegalArgumentException ("a >= b");
      p = alpha;
      q = beta;
      a = aa;
      b = bb;
   }
}