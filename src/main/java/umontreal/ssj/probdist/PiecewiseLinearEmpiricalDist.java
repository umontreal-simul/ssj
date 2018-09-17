/*
 * Class:        PiecewiseLinearEmpiricalDist
 * Description:  piecewise-linear empirical distribution
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

import java.util.Formatter;
import java.util.Locale;
import umontreal.ssj.util.Num;
import umontreal.ssj.util.PrintfFormat;
import java.util.Arrays;
import java.io.IOException;
import java.io.Reader;
import java.io.BufferedReader;

/**
 * Extends the class  @ref ContinuousDistribution for a piecewise-linear
 * approximation of the *empirical* distribution function, based on the
 * observations @f$X_{(1)},…,X_{(n)}@f$ (sorted by increasing order), and
 * defined as follows (e.g., @cite sLAW00a&thinsp; (page 318)). The
 * distribution function starts at @f$X_{(1)}@f$ and climbs linearly by
 * @f$1/(n-1)@f$ between any two successive observations. The density is
 * @f[
 *   f(x) = \frac{1}{(n-1)(X_{(i+1)} - X_{(i)})} \mbox{ for }X_{(i)}\le x < X_{(i+1)}\mbox{ and } i=1,2,…,n-1.
 * @f]
 * The distribution function is
 * @f[
 *   F(x) = \left\{\begin{array}{ll}
 *    0 
 *    & 
 *    \mbox{ for } x < X_{(1)}, 
 *    \\ 
 *    \displaystyle\frac{i-1}{n-1} + \frac{x - X_{(i)}}{(n-1)(X_{(i+1)} - X_{(i)})}
 *    & 
 *   \mbox{ for } X_{(i)} \le x < X_{(i+1)} \mbox{ and } i<n, 
 *    \\ 
 *    1 
 *    & 
 *    \mbox{ for } x \ge X_{(n)}, 
 *   \end{array}\right.
 * @f]
 * whose inverse is
 * @f[
 *   F^{-1}(u) = X_{(i)} + ((n-1)u - i + 1)(X_{(i+1)} - X_{(i)})
 * @f]
 * for @f$(i-1)/(n-1)\le u \le i/(n-1)@f$ and @f$i=1,…,n-1@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class PiecewiseLinearEmpiricalDist extends ContinuousDistribution {
   private double[] sortedObs;
   private double[] diffObs;
   private int n = 0;
   private double sampleMean;
   private double sampleVariance;
   private double sampleStandardDeviation;

/**
 * Constructs a new piecewise-linear distribution using all the observations
 * stored in `obs`. These observations are copied into an internal array and
 * then sorted.
 */
public PiecewiseLinearEmpiricalDist (double[] obs) {
      if (obs.length <= 1)
         throw new IllegalArgumentException ("Two or more observations are needed");
      // sortedObs = obs;
      n = obs.length;
      sortedObs = new double[n];
      System.arraycopy (obs, 0, sortedObs, 0, n);
      init();
   }

   /**
    * Constructs a new empirical distribution using the observations read
    * from the reader `in`. This constructor will read the first `double`
    * of each line in the stream. Any line that does not start with a `+`,
    * `-`, or a decimal digit, is ignored. The file is read until its end.
    * One must be careful about lines starting with a blank. This format
    * is the same as in UNURAN.
    */
   public PiecewiseLinearEmpiricalDist (Reader in) throws IOException {
      BufferedReader inb = new BufferedReader (in);
      double[] data = new double[5];
      String li;
      while ((li = inb.readLine()) != null) {
         // look for the first non-digit character on the read line
         int index = 0;
         while (index < li.length() &&
            (li.charAt (index) == '+' || li.charAt (index) == '-' ||
             li.charAt (index) == 'e' || li.charAt (index) == 'E' ||
             li.charAt (index) == '.' || Character.isDigit (li.charAt (index))))
           ++index; 

         // truncate the line
         li = li.substring (0, index);
         if (!li.equals ("")) {
            try {
               data[n++] = Double.parseDouble (li);
               if (n >= data.length) {
                  double[] newData = new double[2*n];
                  System.arraycopy (data, 0, newData, 0, data.length);
                  data = newData;
               }
            }
            catch (NumberFormatException nfe) {}
         }
      }
      sortedObs = new double[n];
      System.arraycopy (data, 0, sortedObs, 0, n);
      init();
   }


   public double density (double x) {
      // This is implemented via a linear search: very inefficient!!!
      if (x < sortedObs[0] || x >= sortedObs[n-1])
         return 0;
      for (int i = 0; i < (n-1); i++) {
         if (x >= sortedObs[i] && x < sortedObs[i+1])
            return 1.0 / ((n-1)*diffObs[i]);
      }
      throw new IllegalStateException();
   }

   public double cdf (double x) {
      // This is implemented via a linear search: very inefficient!!!
      if (x <= sortedObs[0])
         return 0;
      if (x >= sortedObs[n-1])
         return 1;
      for (int i = 0; i < (n-1); i++) {
         if (x >= sortedObs[i] && x < sortedObs[i+1])
            return i/(n-1.0) + (x - sortedObs[i])/((n-1.0)*diffObs[i]);
      }
      throw new IllegalStateException();
   }

   public double barF (double x) {
      // This is implemented via a linear search: very inefficient!!!
      if (x <= sortedObs[0])
         return 1;
      if (x >= sortedObs[n-1])
         return 0;
      for (int i = 0; i < (n-1); i++) {
         if (x >= sortedObs[i] && x < sortedObs[i+1])
            return (n-1.0-i)/(n-1.0) - (x - sortedObs[i])/((n-1.0)*diffObs[i]);
      }
      throw new IllegalStateException();
   }

   public double inverseF (double u) {
      if (u < 0 || u > 1)
         throw new IllegalArgumentException ("u is not in [0,1]");
      if (u <= 0.0)
         return sortedObs[0];
      if (u >= 1.0)
         return sortedObs[n-1];
      double p = (n - 1)*u;
      int i = (int)Math.floor (p);
      if (i == (n-1))
         return sortedObs[n-1];
      else
         return sortedObs[i] + (p - i)*diffObs[i];
   }

   public double getMean() {
      return sampleMean;
   }

   public double getVariance() {
      return sampleVariance;
   }

   public double getStandardDeviation() {
      return sampleStandardDeviation;
   }

   private void init() {
      Arrays.sort (sortedObs);
      // n = sortedObs.length;
      diffObs = new double[sortedObs.length];
      double sum = 0.0;
      for (int i = 0; i < diffObs.length-1; i++) {
         diffObs[i] = sortedObs[i+1] - sortedObs[i];
         sum += sortedObs[i];
      }
      diffObs[n-1] = 0.0;  // Can be useful in case i=n-1 in inverseF.
      sum += sortedObs[n-1];
      sampleMean = sum / n;
      sum = 0.0;
      for (int i = 0; i < n; i++) {
         double coeff = (sortedObs[i] - sampleMean);
         sum += coeff*coeff;
      }
      sampleVariance = sum / (n-1);
      sampleStandardDeviation = Math.sqrt (sampleVariance);
      supportA = sortedObs[0]*(1.0 - Num.DBL_EPSILON);
      supportB = sortedObs[n-1]*(1.0 + Num.DBL_EPSILON);
   }

/**
 * Returns @f$n@f$, the number of observations.
 */
public int getN() {
      return n;
   }

   /**
    * Returns the value of @f$X_{(i)}@f$.
    */
   public double getObs (int i) {
      return sortedObs[i];
   }

   /**
    * Returns the sample mean of the observations.
    */
   public double getSampleMean() {
      return sampleMean;
   }

   /**
    * Returns the sample variance of the observations.
    */
   public double getSampleVariance() {
      return sampleVariance;
   }

   /**
    * Returns the sample standard deviation of the observations.
    */
   public double getSampleStandardDeviation() {
      return sampleStandardDeviation;
   }

   /**
    * Return a table containing parameters of the current distribution.
    */
   public double[] getParams () {
      double[] retour = new double[n];
      System.arraycopy (sortedObs, 0, retour, 0, n);
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      StringBuilder sb = new StringBuilder();
      Formatter formatter = new Formatter(sb, Locale.US);
      formatter.format(getClass().getSimpleName() + PrintfFormat.NEWLINE);
      for(int i = 0; i<n; i++) {
         formatter.format("%f%n", sortedObs[i]);
      }
      return sb.toString();
   }

}