/*
 * Class:        Pearson6Gen
 * Description:  random variate generators for the Pearson type VI distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements random variate generators for the *Pearson type VI*
 * distribution with shape parameters @f$\alpha_1 > 0@f$ and @f$\alpha_2 >
 * 0@f$, and scale parameter @f$\beta> 0@f$. The density function of this
 * distribution is
 * @anchor REF_randvar_Pearson6Gen_eq_fpearson6
 * @f[
 *   f(x) =\left\{\begin{array}{ll}
 *    \displaystyle\frac{\left(x/{\beta}\right)^{\alpha_1 - 1}}{\beta\mathcal{B}(\alpha_1, \alpha_2)(1 + x/{\beta})^{\alpha_1 + \alpha_2}} 
 *    & 
 *    \quad\mbox{for } x > 0, 
 *    \\ 
 *    0 
 *    & 
 *    \quad\mbox{otherwise,} 
 *   \end{array} \right. \tag{fpearson6}
 * @f]
 * where @f$\mathcal{B}@f$ is the beta function.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class Pearson6Gen extends RandomVariateGen {
   protected double alpha1;
   protected double alpha2;
   protected double beta;

   /**
    * Creates a Pearson6 random variate generator with parameters
    * @f$\alpha_1@f$ = `alpha1`, @f$\alpha_2@f$ = `alpha2` and
    * @f$\beta@f$ = `beta`, using stream `s`.
    */
   public Pearson6Gen (RandomStream s, double alpha1, double alpha2,
                                       double beta) {
      super (s, new Pearson6Dist(alpha1, alpha2, beta));
      setParams (alpha1, alpha2, beta);
   }

   /**
    * Creates a Pearson6 random variate generator with parameters
    * @f$\alpha_1 =@f$ `alpha1`, @f$\alpha_2@f$ = `alpha2` and
    * @f$\beta=1@f$, using stream `s`.
    */
   public Pearson6Gen (RandomStream s, double alpha1, double alpha2) {
      this (s, alpha1, alpha2, 1.0);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public Pearson6Gen (RandomStream s, Pearson6Dist dist) {
      super (s, dist);
      if (dist != null)
         setParams(dist.getAlpha1(), dist.getAlpha2(), dist.getBeta());
   }

   /**
    * Generates a variate from the Pearson VI distribution with shape
    * parameters @f$\alpha_1 > 0@f$ and @f$\alpha_2 > 0@f$, and scale
    * parameter @f$\beta> 0@f$.
    */
   public static double nextDouble (RandomStream s, double alpha1,
                                    double alpha2, double beta) {
      return Pearson6Dist.inverseF (alpha1, alpha2, beta, s.nextDouble());
   }

   /**
    * Returns the @f$\alpha_1@f$ parameter of this object.
    */
   public double getAlpha1() {
      return alpha1;
   }

   /**
    * Returns the @f$\alpha_2@f$ parameter of this object.
    */
   public double getAlpha2() {
      return alpha2;
   }

   /**
    * Returns the @f$\beta@f$ parameter of this object.
    */
   public double getBeta() {
      return beta;
   }

   /**
    * Sets the parameters @f$\alpha_1@f$, @f$\alpha_2@f$ and
    * @f$\beta@f$ of this object.
    */
   public void setParams (double alpha1, double alpha2, double beta) {
      if (alpha1 <= 0.0)
         throw new IllegalArgumentException("alpha1 <= 0");
      if (alpha2 <= 0.0)
         throw new IllegalArgumentException("alpha2 <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException("beta <= 0");
      this.alpha1 = alpha1;
      this.alpha2 = alpha2;
      this.beta = beta;
   }
}