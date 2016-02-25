/*
 * Class:        FBar
 * Description:  Complementary empirical distributions
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
package umontreal.ssj.gof;
import umontreal.ssj.probdist.*;

/**
 * This class is similar to  @ref FDist, except that it provides static
 * methods to compute or approximate the complementary distribution function
 * of @f$X@f$, which we define as @f$\bar{F} (x) = P[X\ge x]@f$, instead of
 * @f$F (x)=P[X\le x]@f$. Note that with our definition of @f$\bar{F}@f$, one
 * has @f$\bar{F} (x) = 1 - F (x)@f$ for continuous distributions and
 * @f$\bar{F} (x) = 1 - F (x-1)@f$ for discrete distributions over the
 * integers.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class FBar {
   private FBar() {}

   private static final double EPSILONSCAN = 1.0E-7;

   private static double scanGlaz (int n, double d, int m) {
      int j, jmoy;
      double temp;
      double jr, jm1r, nr = n;
      int signe;
      double q = 1.0 - d;
      double q4, q3, q2, q1;
      double bin, binMoy;

      jmoy = (int)((n + 1)*d);              // max term of the Binomial
      if (jmoy < m - 1)
         jmoy = m - 1;

      /*---------------------------------------------------------
       * Compute q1: formula (2.5) in Glaz (1989)
       * Compute q2: formula (A.6) in Berman and Eagleson (1985)
       * Compute q3, q4 : Theorem (3.2) in Glaz (1989)
       *---------------------------------------------------------*/

      // compute the probability of term j = jmoy
      q1 = 0.0;
      for (j = 1; j <= jmoy; j++) {
         jr = j;
         q1 += Math.log (nr - jr + 1.0) - Math.log (jr);
      }
      q1 += jmoy*Math.log (d) + (nr - jmoy)*Math.log (q);
      binMoy = Math.exp (q1);
      q1 = binMoy;
      jm1r = jmoy - m + 1;
      if (((jmoy - m + 1) & 1) != 0)
         signe = -1;
      else
         signe = 1;
      q2 = signe*binMoy;
      q3 = signe*binMoy*(2.0 - jm1r*jm1r + jm1r);
      q4 = signe*binMoy*(jm1r + 1.0)*(jm1r + 2.0)*(6.0 + jm1r*jm1r -
         5.0*jm1r);

      // compute the probability of terms j > jmoy
      if (((jmoy - m + 1) & 1) != 0)
         signe = -1;
      else
         signe = 1;

      jm1r = jmoy - m + 1;
      bin = binMoy;
      for (j = jmoy + 1; j <= n; j++) {
         jr = j;
         jm1r += 1.0;
         signe = -signe;
         bin = (bin*(nr - jr + 1.0)*d)/(jr*q);
         if (bin < EPSILONSCAN)
            break;
         q1 += bin;
         q2 += signe*bin;
         q3 += signe*bin*(2.0 - jm1r*jm1r + jm1r);
         q4 += signe*bin*(jm1r + 1.0)*(jm1r + 2.0)*(6.0 + jm1r*jm1r -
            5.0*jm1r);
      }

      q1 = 1.0 - q1;
      q3 /= 2.0;
      q4 /= 12.0;
      if (m == 3) {
        // Problem with this formula; I do not get the same results as Glaz
         q4 = ((nr*(nr - 1.0)*d*d*Math.pow (q, nr - 2.0))/8.0
            + nr*d*2.0*Math.pow (1.0 - 2.0*d, nr - 1.0))
            - 4.0*Math.pow (1.0 - 2.0*d, nr);
         if (d < 1.0/3.0) {
            q4 += nr*d*2.0*Math.pow (1.0 - 3.0*d, nr - 1.0)
                  + 4.0*Math.pow (1.0 - 3.0*d, nr);
         }
      }
      // compute probability: Glaz, equations (3.2) and (3.3)
      q3 = q1 - q2 - q3;
      q4 = q3 - q4;
      //when the approximation is bad, avoid overflow
      temp = Math.log (q3) + (nr - m - 2.0)*Math.log (q4/q3);
      if (temp >= 0.0)
         return 0.0;
      if (temp < -30.0)
         return 1.0;
      q4 = Math.exp (temp);
      return 1.0 - q4;
     }

     //-----------------------------------------------------------------

     private static double scanWNeff (int n, double d, int m) {
      double q = 1.0 - d;
      double temp;
      double bin;
      double sum;
      int j;

      /*--------------------------------------
       * Anderson-Titterington: equation (4)
       *--------------------------------------*/

      // compute the probability of term j = m
      sum = 0.0;
      for (j = 1; j <= m; j++)
         sum += Math.log ((double)(n - j + 1)) - Math.log ((double)j);

      sum += m*Math.log (d) + (n - m)*Math.log (q);
      bin = Math.exp (sum);
      temp = (m/d - n - 1.0)*bin;
      sum = bin;

      // compute the probability of terms j > m
      for (j = m + 1; j <= n; j++) {
         bin *= (n - j + 1)*d/(j*q);
         if (bin < EPSILONSCAN)
            break;
         sum += bin;
      }
      sum = 2.0*sum + temp;
      return sum;
     }

     //----------------------------------------------------------------

     private static double scanAsympt (int n, double d, int m) {
      double kappa;
      double temp;
      double theta;
      double sum;

      /*--------------------------------------------------------------
       * Anderson-Titterington: asymptotic formula after equation (4)
       *--------------------------------------------------------------*/

      theta = Math.sqrt (d/(1.0 - d));
      temp = Math.sqrt ((double)n);
      kappa = m/(d*temp) - temp;
      temp = theta*kappa;
      temp = temp*temp/2.0;
      sum = 2.0*(1.0 - NormalDist.cdf01 (theta*kappa)) +
         (kappa*theta*Math.exp (-temp))/(d*Math.sqrt (2.0*Math.PI));
      return sum;
   }

/**
 * Return @f$P[S_N (d) \ge m]@f$, where @f$S_N (d)@f$ is the scan
 * statistic(see @cite tGLA89a, @cite tGLA01a&thinsp; and  GofStat.scan ),
 * defined as
 * @anchor REF_gof_FBar_eq_scan
 * @f[
 *   S_N (d) = \sup_{0\le y\le1-d} \eta[y, y+d], \tag{scan}
 * @f]
 * where @f$d@f$ is a constant in @f$(0, 1)@f$, @f$\eta[y, y+d]@f$ is the
 * number of observations falling inside the interval @f$[y, y+d]@f$, from a
 * sample of @f$N@f$ i.i.d. @f$U (0,1)@f$ random variables.
 *  One has (see @cite tAND95b&thinsp;),
 * @anchor REF_gof_FBar_DistScan1
 * @anchor REF_gof_FBar_DistScan2
 * @f{align}{
 *    P[S_N (d) \ge m] 
 *    & 
 *   \approx
 *    \left(\frac{m}{d}-N-1\right) b (m) + 2 \sum_{i=m}^N b (i) \tag{DistScan1} 
 *    \\  & 
 *   \approx
 *    2 (1-\Phi(\theta\kappa)) + \theta\kappa\frac{\exp(-\theta^2\kappa^2 /2)}{d \sqrt{2\pi}} \tag{DistScan2}
 * @f}
 * where @f$\Phi@f$ is the standard normal distribution function.
 * @f{align*}{
 *    b (i) 
 *    & 
 *   =
 *    \binom{N}{i} d^i (1-d)^{N-i}, 
 *    \\ 
 *    \theta
 *    & 
 *   =
 *    \sqrt{\frac{d}{1-d}}, 
 *    \\ 
 *    \kappa
 *    & 
 *   =
 *    \frac{m}{d \sqrt{N}} - \sqrt{N}.
 * @f}
 * For @f$d \le1/2@f$, ( {@link REF_gof_FBar_DistScan1
 * DistScan1} ) is exact for @f$m > N/2@f$, but only an approximation
 * otherwise. The approximation ( {@link REF_gof_FBar_DistScan2
 * DistScan2} ) is good when @f$N d^2@f$ is large or when @f$d > 0.3@f$ and
 * @f$N>50@f$. In other cases, this implementation sometimes use the
 * approximation proposed by Glaz @cite tGLA89a&thinsp;. For more
 * information, see @cite tAND95b, @cite tGLA89a, @cite tWAL87a&thinsp;.
 *  The approximation returned by this function is generally good when it is
 * close to 0, but is not very reliable when it exceeds, say, 0.4.  If @f$m
 * \le(N + 1)d@f$, the method returns 1. Else, if @f$Nd \le10@f$, it
 * returns the approximation given by Glaz @cite tGLA89a&thinsp;. If @f$Nd >
 * 10@f$, it computes ( {@link REF_gof_FBar_DistScan2
 * DistScan2} ) or ( {@link REF_gof_FBar_DistScan1 DistScan1} )
 * and returns the result if it does not exceed 0.4, otherwise it computes
 * the approximation from @cite tGLA89a&thinsp;, returns it if it is less
 * than 1.0, and returns 1.0 otherwise.    The relative error can reach 10%
 * when @f$Nd \le10@f$ or when the returned value is less than 0.4. For @f$m
 * > Nd@f$ and @f$Nd > 10@f$, a returned value that exceeds @f$0.4@f$ should
 * be regarded as unreliable. For @f$m = 3@f$, the returned values are
 * totally unreliable. (There may be an error in the original formulae in
 * @cite tGLA89a&thinsp;).  Restrictions: @f$N \ge2@f$ and @f$d
 * \le1/2@f$.<br>
 *  @param n            sample size (@f$\ge2@f$)
 *  @param d            length of the test interval (@f$\in(0,1)@f$)
 *  @param m            scan statistic
 *  @return the complementary distribution function of the statistic
 * evaluated at `m`
 */
public static double scan (int n, double d, int m) {
      double mu;
      double prob;

      if (n < 2)
        throw new IllegalArgumentException ("Calling scan with n < 2");
      if (d <= 0.0 || d >= 1.0)
         throw new IllegalArgumentException ("Calling scan with "+
                    "d outside (0,1)");

      if (m > n)
         return 0.0;
      if (m <= 1)
         return 1.0;
      if (m <= 2) {
         if ((n - 1)*d >= 1.0)
            return 1.0;
         return 1.0 - Math.pow (1.0 - (n - 1)*d, (double)n);
      }
      if (d >= 0.5 && m <= (n + 1)/2.0)
         return 1.0;
      if (d > 0.5)
        return -1.0;              // Error
      // util_Assert (d <= 0.5, "Calling fbar_Scan with d > 1/2");

      mu = n*d;                    // mean of a binomial
      if (m <= mu + d)
         return 1.0;
      if (mu <= 10.0)
         return scanGlaz (n, d, m);
      prob = scanAsympt (n, d, m);
      if ((d >= 0.3 && n >= 50.0) || (n*d*d >= 250.0 && d < 0.3)) {
         if (prob <= 0.4)
            return prob;
      }
      prob = scanWNeff (n, d, m);
      if (prob <= 0.4)
         return prob;
      prob = scanGlaz (n, d, m);
      if (prob > 0.4 && prob <= 1.0)
         return prob;
      return 1.0;
   }
}