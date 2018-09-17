/*
 * Class:        ChiSquareNoncentralDist
 * Description:  Non-central chi-square distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        March 2008
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
import umontreal.ssj.util.*;
import umontreal.ssj.functions.MathFunction;

/**
 * Extends the class  @ref ContinuousDistribution for the *noncentral
 * chi-square* distribution with @f$\nu@f$ degrees of freedom and
 * noncentrality parameter @f$\lambda@f$, where @f$\nu> 0@f$ and
 * @f$\lambda> 0@f$ @cite tJOH95b&thinsp; (page 436). Its density is
 * @anchor REF_probdist_ChiSquareNoncentralDist_eq_nc_fchi2
 * @f[
 *   f(x) =\frac{e^{-(x + \lambda)/2}}{2} \left(\frac{x}{\lambda}\right)^{(\nu-2)/4} I_{\nu/2 - 1}\left(\sqrt{\lambda x}\right) \qquad\mbox{for } x > 0, \tag{nc-fchi2}
 * @f]
 * where @f$I_{\nu}(x)@f$ is the modified Bessel function of the first kind
 * of order @f$\nu@f$ given by
 * @f[
 *   I_{\nu}(z) = \sum_{j=0}^{\infty}\frac{(z/2)^{\nu+ 2j}}{j!\; \Gamma(\nu+ j +1)},
 * @f]
 * where @f$\Gamma(x)@f$ is the gamma function. Notice that this
 * distribution is more general than the *chi-square* distribution since its
 * number of degrees of freedom can be any positive real number. For
 * @f$\lambda= 0@f$ and @f$\nu@f$ a positive integer, we have the ordinary
 * *chi-square* distribution.
 *
 * The cumulative probability function can be written as
 * @anchor REF_probdist_ChiSquareNoncentralDist_eq_nc_cdfchi2
 * @f[
 *   P[X \le x] = \sum_{j=0}^{\infty}\frac{e^{-\lambda/2}(\lambda/2)^j}{j!} P[\chi^2_{\nu+ 2j} \le x], \tag{nc-cdfchi2}
 * @f]
 * where @f$\chi^2_{\nu+ 2j}@f$ is the *central chi-square* distribution
 * with @f$\nu+ 2j@f$ degrees of freedom.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ChiSquareNoncentralDist extends ContinuousDistribution {

   protected static final int PREC = 15;   // num. decimal digits of precision
   protected static final double EPS = Num.TEN_NEG_POW[PREC];

   // To get better precision for very small probabilities
   private static final double PROB_MIN = Double.MIN_VALUE / EPS;
   private static final boolean DETAIL = false;   // For debugging
   private static final double PARLIM = 1000.0;   // Parameters limit

   private double nu;
   private double lambda;
   private PoissonDist pois;
   private int jmin = -1;
   private int jmax = -1;

   private static class Function implements MathFunction {
      protected double nu;
      protected double lambda;
      protected double u;

      public Function (double nu, double lambda, double u) {
         this.nu = nu;
         this.lambda = lambda;
         this.u = u;
      }

      public double evaluate (double x) {
         return u - cdf(nu, lambda, x);
      }
   }

   private static double getXLIM (double nu, double lambda) {
      // for x >= XLIM, density(x) = 0, and cdf(x) = 1
      return (1600.0 + 20.0*(nu + lambda));
   }

   private static int calcJmin (PoissonDist pois) {
      /* Find first lower j such that pois.cdf(j) < EPS. If EPS is chosen
         larger, must change the 7*sig below to find jmin. */
      double lam = pois.getLambda();
      double sig = Math.sqrt(lam);
      int j = (int) (lam - 7.0*sig);
      if (j < 0) return 0;
      while (j > 0 && pois.cdf(j) > EPS)
         j--;
      return j + 1;
   }


   private static int calcJmax (PoissonDist pois) {
      /* Find first upper j such that pois.barF(j) < EPS. If EPS is chosen
         larger, must change the 7*sig below to find jmax. */
      double lam = pois.getLambda();
      double sig = Math.sqrt(lam);
      int j = (int) (lam + 7.0*sig);
      while (pois.barF(j) > EPS)
         j++;
      return j - 1;
   }

   /**
    * Constructs a noncentral chi-square distribution with @f$\nu= @f$
    * `nu` degrees of freedom and noncentrality parameter @f$\lambda= @f$
    * `lambda`.
    */
   public ChiSquareNoncentralDist (double nu, double lambda) {
      setParams (nu, lambda);
   }


   public double density (double x) {
      return density (nu, lambda, x);
   }

   public double cdf (double x) {
      if (x <= 0.0)
         return 0.0;
      if (x >= getXLIM (nu, lambda))
         return 1.0;

      if (nu >= PARLIM || lambda >= PARLIM)
         return cdfPenev (false, nu, lambda, x);

      return cdfExact (pois, jmax, nu, lambda, x);
   }

   public double barF (double x) {
      if (x <= 0.0)
         return 1.0;
      if (x >= getXLIM (nu, lambda))
         return 0.0;

      if (nu >= PARLIM || lambda >= PARLIM)
         return cdfPenev (true, nu, lambda, x);

      return barFExact (pois, jmin, nu, lambda, x);
   }

   public double inverseF (double u) {
      return inverseF (nu, lambda, u);
   }

   public double getMean() {
      return ChiSquareNoncentralDist.getMean (nu, lambda);
   }

   public double getVariance() {
      return ChiSquareNoncentralDist.getVariance (nu, lambda);
   }

   public double getStandardDeviation() {
      return ChiSquareNoncentralDist.getStandardDeviation (nu, lambda);
   }

/**
 * Computes the density function (
 * {@link REF_probdist_ChiSquareNoncentralDist_eq_nc_fchi2
 * nc-fchi2} ) for a *noncentral chi-square* distribution with @f$\nu= @f$
 * `nu` degrees of freedom and parameter @f$\lambda= @f$ `lambda`.
 */
public static double density (double nu, double lambda, double x) {
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (x <= 0.0)
         return 0.0;
      if (lambda == 0.0)
         return GammaDist.density (nu/2.0, 0.5, x);
      if (x >= getXLIM (nu, lambda))
         return 0.0;

      final double a = 0.5*nu - 1.0;      // a = order of Bessel I_a(x)
      double z = lambda*x;

      // The series term is maximal for j = JMAX
      final long JMAX = (long) (0.5*z/(a + Math.sqrt(a*a + z)));

      // calculer le terme pour j = JMAX --> termMax
      double termMax = -0.5*(x + lambda) + 0.25*(nu - 2)* Math.log (x/lambda) +
         0.5*(a + 2*JMAX)*Math.log(0.25*x * lambda) - Num.lnFactorial(JMAX) -
         Num.lnGamma(a + JMAX + 1) - Num.LN2;

      termMax = Math.exp(termMax);

      // Sum non negligible terms on each side of the max term
      double sum = termMax;
      double term = termMax;
      long j = JMAX;
      double y = 4.0/z;

      while (j > 0 && term > sum*EPS) {
         term *= y*j*(j + a);
         sum += term;
         j--;
      }

      term = termMax;
      j = JMAX + 1;
      y = z/4.0;
      while (term > sum*EPS) {
         term *= y/(j*(j + a));
         sum += term;
         j++;
      }
      return sum;
   }

   /**
    * Computes the noncentral chi-square distribution function (
    * {@link REF_probdist_ChiSquareNoncentralDist_eq_nc_cdfchi2
    * nc-cdfchi2} ) with @f$\nu= @f$ `nu` degrees of freedom and
    * parameter @f$\lambda= @f$ `lambda`.  For @f$\lambda< 1000@f$, we
    * use a reasonnably efficient method of summing the cumulative
    * function (
    * {@link REF_probdist_ChiSquareNoncentralDist_eq_nc_cdfchi2
    * nc-cdfchi2} ) by using recurrence relations to compute both the
    * Poisson and the chi-square probabilities in terms of neighboring
    * probability terms. For @f$\lambda\ge1000@f$, we use the normal
    * approximation given in @cite tPEN00a&thinsp;, except very near the
    * mean where we use the central chi square approximation from
    * @cite tPEA59a&thinsp;. This function returns at least 6 decimal
    * degits of precision nearly everywhere.
    */
   public static double cdf (double nu, double lambda, double x) {
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (x <= 0.0)
         return 0.0;
      if (lambda == 0.0)
         return GammaDist.cdf (nu/2.0, 0.5, PREC, x);
      if (x >= getXLIM (nu, lambda))
         return 1.0;

      if (nu >= PARLIM || lambda >= PARLIM)
         return cdfPenev (false, nu, lambda, x);

      PoissonDist pois = new PoissonDist (lambda / 2.0);
      int j = calcJmax (pois);
      return cdfExact (pois, j, nu, lambda, x);
   }


   private static double cdfExact (PoissonDist pois, int jmax, double nu,
         double lambda, double x) {
      final int JMED = (int) (lambda / 2.0);

      int j = jmax;
      double chicdf = GammaDist.cdf (j + 0.5*nu, 0.5, PREC, x);
      if (chicdf >= 1.0) return 1.0;

      double chiterm = x/(j + 0.5*nu) * GammaDist.density(j + 0.5*nu, 0.5, x);
      double prob = pois.prob(j);
      double sum = prob * chicdf;

      if (DETAIL) {
         System.out.println (PrintfFormat.NEWLINE +
             "----------------------------------- cdf efficace");
         System.out.print ("   j            chi term                chi CDF");
         System.out.println ("                   Poisson p            somme" +
             PrintfFormat.NEWLINE);
         System.out.println (PrintfFormat.d (5, j) + "   " +
                             PrintfFormat.g (20, 10, chiterm) + "   " +
                             PrintfFormat.g (22, 15, chicdf) + "   " +
                             PrintfFormat.g (20, 10, prob) + "   " +
                             PrintfFormat.g (22, 15, sum));
      }

      --j;
      while (j >= 0 && (prob > sum*EPS || j >= JMED)) {
         if (chicdf <= PROB_MIN) {
             chicdf = GammaDist.cdf (j + nu/2.0, 0.5, PREC, x);
             chiterm = x / (j + 0.5*nu) *GammaDist.density(j + nu/2.0, 0.5, x);
         } else {
            chiterm *= (2 + 2*j + nu) / x;
            chicdf += chiterm;
         }
         if (chicdf >= 1.0 - EPS)
            return sum + pois.cdf(j);
         prob = pois.prob(j);
         sum += prob * chicdf;
         if (DETAIL) {
            System.out.println (PrintfFormat.d (5, j) + "   " +
                                PrintfFormat.g (20, 10, chiterm) + "   " +
                                PrintfFormat.g (22, 15, chicdf) + "   " +
                                PrintfFormat.g (20, 10, prob) + "   " +
                                PrintfFormat.g (22, 15, sum));
         }
         --j;
      }
      return sum;
   }

/*******************
   public static double cdf4 (double nu, double lambda, double x)
   {
      // Version naive et lente
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (x <= 0.0)
         return 0.0;
      if (lambda == 0.0)
         return GammaDist.cdf (nu/2.0, 0.5, PREC, x);
      PoissonDist pois = new PoissonDist (lambda / 2.0);
      int jmed = (int) (lambda / 2.0);
      int j = 0;
      double sum = 0.0;
      double prob = 1.0e100;

      if (DETAIL) {
         System.out.println (PrintfFormat.NEWLINE +
             "-------------------------------------- cdf naif");
         System.out.print ("   j             Poisson p               chi CDF");
         System.out.println ("            somme" + PrintfFormat.NEWLINE);
      }

      while (j <= jmed || prob > sum*EPS)
      {
         double chi2cdf = GammaDist.cdf (j + nu/2.0, 0.5, PREC, x);
         prob = pois.prob(j);
         sum += prob * chi2cdf;
         if (DETAIL)
            System.out.println (PrintfFormat.d (5, j) + "   " +
                                PrintfFormat.g (20, 10, prob) + "   " +
                                PrintfFormat.g (20, 10, chi2cdf) + "   " +
                                PrintfFormat.g (20, 10, sum));
         ++j;
      }
      if (DETAIL)
         System.out.println ("");
      return sum;
   }
********************/

   private static double cdfPearson (boolean bar, double nu, double lambda,
                                     double x) {
      // Pearson's central chi square approximation from @cite tPEA59a
      // If bar is true, then returns u = barF; else returns u = cdf
      double t2 = (nu + 2.0*lambda);
      double t3 = (nu + 3.0*lambda);
      double lib = t2*t2*t2/(t3*t3);
      double z = x + lambda*lambda/t3;
      if (bar)
         return GammaDist.barF(lib/2.0, 0.5, PREC, t2*z/t3);
      else
         return GammaDist.cdf(lib/2.0, 0.5, PREC, t2*z/t3);
   }

   private static double penevH (double y) {
      // The h function in Penev-Raykov paper
      if (0.0 == y)
         return 0.0;
      if (1.0 == y)
         return 0.5;
      return ((1.0 - y)*Math.log1p(-y) + y - 0.5*y*y) / (y*y);
   }


   private static double penevB (double mu, double s, double h1) {
      // The B function in Penev-Raykov paper
      final double f = 1.0 + 2.0 * mu * s;
      final double g = 1.0 + 3.0 * mu * s;
      final double c = (f - 2.0 * h1 - s *f) / (f - 2.0 * h1);
      double B = -1.5 * (1.0 + 4.0 * mu * s) / (f * f);
      B += 5.0 * g * g / (3.0 * f * f * f) + 2.0 * g / ((s - 1.0) * f * f);
      B += 3.0 * c / ((s - 1.0) * (s - 1.0) * f);
      B -= c * c * (1.0 + 2.0 * penevH (c)) / (2.0 * (s - 1.0) * (s - 1.0) * f);
      return B;
   }

   private static double getPenevLim (double nu, double lam) {
      if (lam >= 100000.0)
         return 5.0;
      if (nu >= 20000.0)
         return 3.0;
      return 2.0;
   }

   private static double cdfPenev (boolean bar, double nu, double lambda,
                                   double x) {
      /* Penev-Raykov normal approximation from @cite tPEN00a
         If bar is true, then returns u = barF; else returns u = cdf.
         This approximation is very good except very near the mean where it
         is bad or give a NaN. Thus, near the mean, we use the Pearson approx.
      */
      double lim = getPenevLim(nu, lambda);
      double mean = nu + lambda;
      if (x >= mean - lim && x <= mean + lim)
         return cdfPearson (bar, nu, lambda, x);

      final double mu = lambda / nu;
      final double s = (Math.sqrt(1.0 + 4.0 * x * mu / nu) - 1.0) / (2.0 * mu);
      if (s == 1.0)
         return 0.5;
      final double h1 = penevH (1.0 - s);
      final double B = penevB (mu, s, h1);
      final double f = 1.0 + 2.0 * mu * s;
      double z = nu * (s - 1.0) * (s - 1.0) * (0.5 / s + mu - h1 / s) -
                 Math.log(1.0 / s - 2.0 * h1 / (s * f)) + 2.0 * B / nu;
      // If z is a NaN, then z != z
      if (z < 0.0 || z != z)
          return cdfPearson (bar, nu, lambda, x);

      z = Math.sqrt(z);
      if (s < 1)
         z = -z;
      if (bar)
         return NormalDist.barF01(z);
      else
         return NormalDist.cdf01(z);
   }

   /**
    * Computes the complementary noncentral chi-square distribution
    * function with @f$\nu= @f$ `nu` degrees of freedom and parameter
    * @f$\lambda= @f$ `lambda`.
    */
   public static double barF (double nu, double lambda, double x) {
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (x <= 0.0)
         return 1.0;
      if (lambda == 0.0)
         return GammaDist.barF (nu/2.0, 0.5, PREC, x);
      if (x >= getXLIM (nu, lambda))
         return 0.0;

      if (nu >= PARLIM || lambda >= PARLIM)
         return cdfPenev (true, nu, lambda, x);

      PoissonDist pois = new PoissonDist (lambda / 2.0);
      int j = calcJmin (pois);
      return barFExact (pois, j, nu, lambda, x);
   }

/***********************
   public static double bar4 (double nu, double lambda, double x)
   {
      // Version naive et lente
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (x <= 0.0)
         return 1.0;
      if (lambda == 0.0)
         return GammaDist.barF (nu/2.0, 0.5, PREC, x);
      PoissonDist pois = new PoissonDist (lambda / 2.0);
      int jmed = (int) (lambda / 2.0);
      int j = 0;
      double sum = 0.0;
      double pdfterm = 0.0;
      double prob = 1.0e100;

      if (DETAIL) {
         System.out.println (PrintfFormat.NEWLINE +
             "-------------------------------------- barF naif");
         System.out.print ("   j             Poisson p               chi term");
         System.out.println ("               chi barF            somme" +
             PrintfFormat.NEWLINE);
      }

      while (j <= jmed || prob > sum*EPS)
      {
         double chi2cdf = GammaDist.barF (j + nu/2.0, 0.5, PREC, x);
         prob = pois.prob(j);
         sum += prob * chi2cdf;
         if (DETAIL)
            System.out.println (PrintfFormat.d (5, j) + "   " +
                                PrintfFormat.g (20, 12, prob) + "   " +
                       PrintfFormat.g (20, 12, chi2cdf - pdfterm) + "   " +
                                PrintfFormat.g (20, 12, chi2cdf) + "   " +
                                PrintfFormat.g (20, 12, sum));
         ++j;
         pdfterm = chi2cdf;
      }
      if (DETAIL)
         System.out.println ("");
      return sum;
   }
********************/

   private static double barFExact (PoissonDist pois, int jmin, double nu,
         double lambda, double x) {
      final int JMED = (int) (lambda / 2.0);

      int j = jmin;
      double chibar = GammaDist.barF (j + 0.5*nu, 0.5, PREC, x);
      if (chibar >= 1.0) return 1.0;

      double prob = pois.prob(j);
      double sum = prob * chibar;
      double chiterm = 2.0* GammaDist.density(j + 0.5*nu, 0.5, x);

      if (DETAIL) {
         System.out.println (PrintfFormat.NEWLINE +
              "--------------------------------- barF efficace");
         System.out.print ("   j                  Poisson p            chi term");
         System.out.println ("                chi barF                 somme" +
              PrintfFormat.NEWLINE);
         System.out.println (PrintfFormat.d (5, j) + "   " +
                             PrintfFormat.g (20, 12, prob) + "   " +
                             PrintfFormat.g (20, 12, chiterm) + "   " +
                             PrintfFormat.g (22, 12, chibar) + "   " +
                             PrintfFormat.g (22, 12, sum));
      }

      ++j;
      while (prob > sum*EPS || j <= JMED)
      {
         if (chibar <= PROB_MIN) {
             chibar = GammaDist.barF (j + nu/2.0, 0.5, PREC, x);
             chiterm = 2.0 *GammaDist.density(j + 0.5*nu, 0.5, x);
         } else {
            chiterm *= x/(2*j - 2 + nu);
            chibar += chiterm;
         }
         if (chibar >= 1.0)
            return sum + pois.barF(j);
         prob = pois.prob(j);
         sum += prob * chibar;

         if (DETAIL) {
            System.out.println (PrintfFormat.d (5, j) + "   " +
                                PrintfFormat.g (20, 12, prob) + "   " +
                                PrintfFormat.g (20, 12, chiterm) + "   " +
                                PrintfFormat.g (22, 12, chibar) + "   " +
                                PrintfFormat.g (22, 12, sum));
         }
         ++j;
      }
      return sum;
   }

   /**
    * Computes the inverse of the noncentral chi-square distribution with
    * @f$\nu= @f$ `nu` degrees of freedom and parameter @f$\lambda= @f$
    * `lambda`.
    */
   public static double inverseF (double nu, double lambda, double u) {
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;
      if (u <= 0.0)
         return 0.0;
      if (lambda == 0.0)
         return GammaDist.inverseF (nu/2.0, 0.5, PREC, u);

      double x = inverse9 (nu, lambda, u);
      double v = cdf (nu, lambda, x);

      Function f = new Function (nu, lambda, u);
      if (v >= u)
         x = RootFinder.brentDekker (0.0, x, f, 1.0e-10);
      else {
         v = getXLIM (nu, lambda);
         x = RootFinder.brentDekker (x, v, f, 1.0e-10);
      }
      return x;
   }

   private static double inverse9 (double nu, double lambda, double u) {
       // Une approximation normale @cite tKRI06a
       double z = NormalDist.inverseF01 (u);
       double a = nu + lambda;
       double b = (nu + 2.0*lambda)/(a*a);
       double t = z*Math.sqrt(2.0*b/9.0) - 2.0*b/9.0 + 1.0;
       return a*t*t*t;
   }

   /**
    * Computes and returns the mean @f$E[X] = \nu+ \lambda@f$ of the
    * noncentral chi-square distribution with parameters @f$\nu= @f$ `nu`
    * and @f$\lambda= @f$ `lambda`.
    *  @return the mean of the Noncentral noncentral chi-square
    * distribution
    */
   public static double getMean (double nu, double lambda) {
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      return  nu + lambda;
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] = 2(\nu+
    * 2\lambda)@f$ of the noncentral chi-square distribution with
    * parameters @f$\nu=@f$ `nu` and @f$\lambda= @f$ `lambda`.
    *  @return the variance of the noncentral chi-square distribution
    */
   public static double getVariance (double nu, double lambda) {
      if (nu <= 0.)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      return  2.0 * (nu + 2.0*lambda);
   }

   /**
    * Computes and returns the standard deviation of the noncentral
    * chi-square distribution with parameters @f$\nu=@f$ `nu` and
    * @f$\lambda= @f$ `lambda`.
    *  @return the standard deviation of the noncentral chi-square
    * distribution
    */
   public static double getStandardDeviation (double nu, double lambda) {
      return Math.sqrt(getVariance (nu, lambda));
   }

   /**
    * Returns the parameter @f$\nu@f$ of this object.
    */
   public double getNu() {
      return nu;
   }

   /**
    * Returns the parameter @f$\lambda@f$ of this object.
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Sets the parameters @f$\nu=@f$ `nu` and @f$\lambda= @f$ `lambda`
    * of this object.
    */
   public void setParams (double nu, double lambda) {
      if (nu <= 0.0)
         throw new IllegalArgumentException ("nu <= 0");
      if (lambda < 0.0)
         throw new IllegalArgumentException ("lambda < 0");
      if (lambda == 0.0)
         throw new IllegalArgumentException ("lambda = 0");
      this.nu = nu;
      this.lambda = lambda;
      supportA = 0.0;
      if (nu >= PARLIM || lambda >= PARLIM)
         return;
      pois = new PoissonDist (lambda / 2.0);
      jmax = calcJmax (pois);
      jmin = calcJmin (pois);
   }

   /**
    * Returns a table containing the parameters of the current
    * distribution.
    */
   public double[] getParams () {
      double[] retour = {nu, lambda};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + ":   nu = " + nu + ",   lambda = " + lambda;
   }

}