/*
 * Class:        UniformIntDist
 * Description:  discrete uniform distribution
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

/**
 * Extends the class  @ref DiscreteDistributionInt for the *discrete uniform*
 * distribution over the range @f$[i,j]@f$. Its mass function is given by
 * @anchor REF_probdist_UniformIntDist_eq_fmassuniformint
 * @f[
 *   p(x) = \frac{1}{j - i + 1} \qquad\mbox{ for } x = i, i + 1, …, j \tag{fmassuniformint}
 * @f]
 * and 0 elsewhere. The distribution function is
 * @anchor REF_probdist_UniformIntDist_eq_cdfuniformint
 * @f[
 *   F(x) = \left\{\begin{array}{ll}
 *    0, 
 *    & 
 *    \mbox{ for } x < i
 *    \\ 
 *    \displaystyle\frac{\lfloor x\rfloor-i+1}{j-i+1}, 
 *    & 
 *    \mbox{ for } i\le x < j
 *    \\ 
 *    1, 
 *    & 
 *    \mbox{ for } x \ge j. 
 *   \end{array}\right. \tag{cdfuniformint}
 * @f]
 * and its inverse is
 * @anchor REF_probdist_UniformIntDist_eq_invuniformint
 * @f[
 *   F^{-1}(u) = i + \lfloor(j - i + 1)u\rfloor\qquad\mbox{for }0 \le u \le1. \tag{invuniformint}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discrete
 */
public class UniformIntDist extends DiscreteDistributionInt {
   protected int i;
   protected int j;

   /**
    * Constructs a discrete uniform distribution over the interval
    * @f$[i,j]@f$.
    */
   public UniformIntDist (int i, int j) {
      setParams (i, j);
   }


   public double prob (int x) {
      return prob (i, j, x);
   }

   public double cdf (int x) {
      return cdf (i, j, x);
   }

   public double barF (int x) {
      return barF (i, j, x);
   }

   public int inverseFInt (double u) {
      return inverseF (i, j, u);
   }

   public double getMean() {
      return getMean (i, j);
   }

   public double getVariance() {
      return getVariance (i, j);
   }

   public double getStandardDeviation() {
      return getStandardDeviation (i, j);
   }

/**
 * Computes the discrete uniform probability @f$p(x)@f$ defined in (
 * {@link REF_probdist_UniformIntDist_eq_fmassuniformint
 * fmassuniformint} ).
 */
public static double prob (int i, int j, int x) {
      if (j < i)
         throw new IllegalArgumentException ("j < i");
      if (x < i || x > j)
         return 0.0;

      return (1.0 / (j - i + 1.0));
   }

   /**
    * Computes the discrete uniform distribution function defined in (
    * {@link REF_probdist_UniformIntDist_eq_cdfuniformint
    * cdfuniformint} ).
    */
   public static double cdf (int i, int j, int x) {
      if (j < i)
         throw new IllegalArgumentException ("j < i");
      if (x < i)
         return 0.0;
      if (x >= j)
         return 1.0;

      return ((x - i + 1) / (j - i + 1.0));
   }

   /**
    * Computes the discrete uniform complementary distribution function
    * @f$\bar{F}(x)@f$. *WARNING:* The complementary distribution function
    * is defined as @f$\bar{F}(x) = P[X \ge x]@f$.
    */
   public static double barF (int i, int j, int x) {
      if (j < i)
        throw new IllegalArgumentException ("j < i");
      if (x <= i)
         return 1.0;
      if (x > j)
         return 0.0;

      return ((j - x + 1.0) / (j - i + 1.0));
   }

   /**
    * Computes the inverse of the discrete uniform distribution function (
    * {@link REF_probdist_UniformIntDist_eq_invuniformint
    * invuniformint} ).
    */
   public static int inverseF (int i, int j, double u) {
      if (j < i)
        throw new IllegalArgumentException ("j < i");

       if (u > 1.0 || u < 0.0)
           throw new IllegalArgumentException ("u not in [0, 1]");

       if (u <= 0.0)
           return i;
       if (u >= 1.0)
           return j;

       return i + (int) (u * (j - i + 1.0));
   }

   /**
    * Estimates the parameters @f$(i, j)@f$ of the uniform distribution
    * over integers using the maximum likelihood method, from the @f$n@f$
    * observations @f$x[k]@f$, @f$k = 0, 1, …, n-1@f$. The estimates are
    * returned in a two-element array, in regular order: [@f$i@f$,
    * @f$j@f$].  The maximum likelihood estimators are the values
    * @f$(\hat{\imath}@f$, @f$\hat{\jmath})@f$ that satisfy the
    * equations
    * @f{align*}{
    *    \hat{\imath} 
    *    & 
    *    = 
    *    \mbox{min} \{x_k\}
    *    \\ 
    *   \hat{\jmath} 
    *    & 
    *    = 
    *    \mbox{max} \{x_k\}
    * @f}
    * where @f$\bar{x}_n@f$ is the average of @f$x[0],…,x[n-1]@f$.
    *  @param x            the list of observations used to evaluate
    *                      parameters
    *  @param n            the number of observations used to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\imath}@f$,
    * @f$\hat{\jmath}@f$]
    */
   public static double[] getMLE (int[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double parameters[] = new double[2];
      parameters[0] = (double) Integer.MAX_VALUE;
      parameters[1] = (double) Integer.MIN_VALUE;
      for (int i = 0; i < n; i++) {
         if ((double) x[i] < parameters[0])
            parameters[0] = (double) x[i];
         if ((double) x[i] > parameters[1])
            parameters[1] = (double) x[i];
      }
      return parameters;
   }

   /**
    * Creates a new instance of a discrete uniform distribution over
    * integers with parameters @f$i@f$ and @f$j@f$ estimated using the
    * maximum likelihood method based on the @f$n@f$ observations
    * @f$x[k]@f$, @f$k = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static UniformIntDist getInstanceFromMLE (int[] x, int n) {

      double parameters[] = getMLE (x, n);

      return new UniformIntDist ((int) parameters[0], (int) parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = (i + j)/2@f$ of the discrete
    * uniform distribution.
    *  @return the mean of the discrete uniform distribution
    */
   public static double getMean (int i, int j) {
      if (j < i)
        throw new IllegalArgumentException ("j < i");

      return ((i + j) / 2.0);
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = [(j - i +
    * 1)^2 - 1]/{12}@f$ of the discrete uniform distribution.
    *  @return the variance of the discrete uniform distribution
    */
   public static double getVariance (int i, int j) {
      if (j < i)
        throw new IllegalArgumentException ("j < i");

      return (((j - i + 1.0) * (j - i + 1.0) - 1.0) / 12.0);
   }

   /**
    * Computes and returns the standard deviation of the discrete uniform
    * distribution.
    *  @return the standard deviation of the discrete uniform distribution
    */
   public static double getStandardDeviation (int i, int j) {
      return Math.sqrt (UniformIntDist.getVariance (i, j));
   }

   /**
    * Returns the parameter @f$i@f$.
    */
   public int getI() {
      return i;
   }

   /**
    * Returns the parameter @f$j@f$.
    */
   public int getJ() {
      return j;
   }

   /**
    * Sets the parameters @f$i@f$ and @f$j@f$ for this object.
    */
   public void setParams (int i, int j) {
      if (j < i)
        throw new IllegalArgumentException ("j < i");

      supportA = this.i = i;
      supportB = this.j = j;
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in regular order: [@f$i@f$,
    * @f$j@f$].
    */
   public double[] getParams () {
      double[] retour = {i, j};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : i = " + i + ", j = " + j;
   }

}