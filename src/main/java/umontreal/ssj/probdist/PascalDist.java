/*
 * Class:        PascalDist
 * Description:  Pascal distribution
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
package umontreal.ssj.probdist;
import umontreal.ssj.util.RootFinder;
import umontreal.ssj.functions.MathFunction;

/**
 * The *Pascal* distribution is a special case of the *negative binomial*
 * distribution @cite sLAW00a&thinsp; (page 324) with parameters @f$n@f$ and
 * @f$p@f$, where @f$n@f$ is a positive integer and @f$0\le p\le1@f$. Its
 * mass function is
 * @anchor REF_probdist_PascalDist_eq_fmass_pascal
 * @f[
 *   p(x) = \binom{n + x - 1}{x} p^n (1 - p)^x, \qquad\mbox{for } x = 0, 1, 2, …\tag{fmass-pascal}
 * @f]
 * This @f$p(x)@f$ can be interpreted as the probability of having @f$x@f$
 * failures before the @f$n@f$th success in a sequence of independent
 * Bernoulli trials with probability of success @f$p@f$. For @f$n=1@f$, this
 * gives the *geometric* distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class PascalDist extends NegativeBinomialDist {
   private static final double EPSI = 1.0E-10;

   private static class Function implements MathFunction {
      protected int m;
      protected int max;
      protected double mean;
      protected int[] Fj;

      public Function (int m, int max, double mean, int[] Fj) {
         this.m = m;
         this.max = max;
         this.mean = mean;
         this.Fj = new int[Fj.length];
         System.arraycopy(Fj, 0, this.Fj, 0, Fj.length);
      }

      public double evaluate (double p) {
         double sum = 0.0;
         double s = (p * mean) / (1.0 - p);

         for (int j = 0; j < max; j++)
            sum += Fj[j] / (s + (double) j);

         return sum + m * Math.log (p);
      }

      public double evaluateN (int n, double p) {
         double sum = 0.0;

         for (int j = 0; j < max; j++)
            sum += Fj[j] / (n + j);

         return sum + m * Math.log (p);
      }
   }

   /**
    * Creates an object that contains the probability terms (
    * {@link REF_probdist_PascalDist_eq_fmass_pascal
    * fmass-pascal} ) and the distribution function for the Pascal
    * distribution with parameter @f$n@f$ and @f$p@f$.
    */
   public PascalDist (int n, double p) {
      setParams (n, p);
   }

   /**
    * Estimates the parameter @f$(n, p)@f$ of the Pascal distribution
    * using the maximum likelihood method, from the @f$m@f$ observations
    * @f$x[i]@f$, @f$i = 0, 1, …, m-1@f$. The estimates are returned in a
    * two-element array, in regular order: [@f$n@f$, @f$p@f$].  The
    * maximum likelihood estimators are the values @f$(\hat{n}@f$,
    * @f$\hat{p})@f$ that satisfy the equations
    * @f{align*}{
    *    \frac{\hat{n}(1 - \hat{p})}{\hat{p}} 
    *    & 
    *    = 
    *    \bar{x}_m
    *    \\ 
    *   \ln(1 + \hat{p}) 
    *    & 
    *    = 
    *    \sum_{j=1}^{\infty} \frac{F_j}{(\hat{n} + j - 1)}
    * @f}
    * where @f$\bar{x}_m@f$ is the average of @f$x[0],…,x[m-1]@f$, and
    * @f$F_j = \sum_{i=j}^{\infty} f_i@f$ = proportion of @f$x@f$’s
    * which are greater than or equal to @f$j@f$ @cite tJOH69a&thinsp;
    * (page 132).
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param m            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{n}@f$, @f$\hat{p}@f$]
    */
   public static double[] getMLE (int[] x, int m) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");

      double sum = 0.0;
      int max = Integer.MIN_VALUE;
      for (int i = 0; i < m; i++)
      {
         sum += x[i];
         if (x[i] > max)
            max = x[i];
      }
      double mean = (double) sum / (double) m;

      double var = 0.0;
      for (int i = 0; i < m; i++)
         var += (x[i] - mean) * (x[i] - mean);
      var /= (double) m;

      if (mean >= var)
           throw new UnsupportedOperationException("mean >= variance");

      int[] Fj = new int[max];
      for (int j = 0; j < max; j++) {
         int prop = 0;
         for (int i = 0; i < m; i++)
            if (x[i] > j)
               prop++;

         Fj[j] = prop;
      }

      double[] parameters = new double[2];
      Function f = new Function (m, max, mean, Fj);

      parameters[1] = RootFinder.brentDekker (EPSI, 1 - EPSI, f, 1e-5);
      if (parameters[1] >= 1.0)
         parameters[1] = 1.0 - 1e-15;

      parameters[0] = Math.round ((parameters[1] * mean) / (1.0 - parameters[1]));
      if (parameters[0] == 0)
          parameters[0] = 1;

      return parameters;
   }

   /**
    * Creates a new instance of a Pascal distribution with parameters
    * @f$n@f$ and @f$p@f$ estimated using the maximum likelihood method
    * based on the @f$m@f$ observations @f$x[i]@f$, @f$i = 0, 1, …,
    * m-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param m            the number of observations to use to evaluate
    *                      parameters
    */
   public static PascalDist getInstanceFromMLE (int[] x, int m) {
      double parameters[] = getMLE (x, m);
      return new PascalDist ((int) parameters[0], parameters[1]);
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public int getN1() {
      return (int) (n + 0.5);
   }

   /**
    * Sets the parameter @f$n@f$ and @f$p@f$ of this object.
    */
   public void setParams (int n, double p) {
      super.setParams ((double) n, p);
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : n = " + getN1() + ", p = " + getP();
   }

}