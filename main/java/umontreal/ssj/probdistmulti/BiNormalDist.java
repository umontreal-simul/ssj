/*
 * Class:        BiNormalDist
 * Description:  bivariate normal distribution
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
package umontreal.ssj.probdistmulti;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.NormalDistQuick;

/**
 * Extends the class  @ref ContinuousDistribution2Dim for the *bivariate
 * normal* distribution @cite tJOH72a&thinsp; (page 84). It has means @f$E[X]
 * =\mu_1@f$, @f$E[Y] =\mu_2@f$, and variances var@f$[X] = \sigma_1^2@f$,
 * var@f$[Y] = \sigma_2^2@f$ such that @f$\sigma_1 > 0@f$ and @f$\sigma_2
 * > 0@f$. The correlation between @f$X@f$ and @f$Y@f$ is @f$\rho@f$. Its
 * density function is
 * @anchor REF_probdistmulti_BiNormalDist_eq_f1binormal
 * @f[
 *   f (x, y) = \frac{1}{2\pi\sigma_1\sigma_2\sqrt{1-\rho^2}}e^{-T}  \tag{f1binormal}
 * @f]
 * @f[
 *   T = \frac{1}{2(1-\rho^2)}\left[\left(\frac{x-\mu_1}{\sigma_1}\right)^2 - 2\rho\left(\frac{x-\mu_1}{\sigma_1}\right) \left(\frac{y-\mu_2}{\sigma_2}\right) + \left(\frac{y-\mu_2}{\sigma_2}\right)^2\right]
 * @f]
 * and the corresponding distribution function is (the `cdf` method)
 * @anchor REF_probdistmulti_BiNormalDist_eq_cdf1binormal
 * @f[
 *   \Phi(\mu_1, \sigma_1, x, \mu_2, \sigma_2, y, \rho) = \frac{1}{2\pi\sigma_1\sigma_2\sqrt{1-\rho^2}} \int_{-\infty}^x dx \int_{-\infty}^y dy  e^{-T}. \tag{cdf1binormal}
 * @f]
 * We also define the upper distribution function (the `barF` method) as
 * @anchor REF_probdistmulti_BiNormalDist_eq_cdf3binormal
 * @f[
 *   \overline{\Phi}(\mu_1, \sigma_1, x, \mu_2, \sigma_2, y, \rho) = \frac{1}{2\pi\sigma_1\sigma_2\sqrt{1-\rho^2}} \int^{\infty}_x dx \int^{\infty}_y dy  e^{-T}. \tag{cdf3binormal}
 * @f]
 * When @f$\mu_1=\mu_2=0@f$ and @f$\sigma_1=\sigma_2=1@f$, we have the
 * *standard binormal* distribution, with corresponding distribution function
 * @anchor REF_probdistmulti_BiNormalDist_eq_cdf2binormal
 * @f[
 *   \Phi(x, y, \rho) = \frac{1}{2\pi\sqrt{1-\rho^2}} \int_{-\infty}^x dx \int_{-\infty}^y dy  e^{-S} \tag{cdf2binormal}
 * @f]
 * @f[
 *   S = \frac{x^2 - 2\rho x y + y^2}{2(1-\rho^2)}.
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_continuous2dim
 */
public class BiNormalDist extends ContinuousDistribution2Dim {
   protected int ndigit;        // Number of decimal digits of accuracy
   protected double mu1, mu2;
   protected double sigma1, sigma2;
   protected double rho;
   protected double racRho;        // sqrt(1 - rho^2)
   protected double detS;          // 2*PI*sigma1*sigma2*sqrt(1 - rho^2)
   protected static final double RHO_SMALL = 1.0e-8; // neglect small rhos

   private static final double Z[] = { 
      0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992 };

   private static final double W[] = {
      0.018854042, 0.038088059, 0.0452707394, 0.038088059, 0.018854042 };

   private static final double AGauss[] = {
      -0.72657601, 0.71070688, -0.142248368, 0.127414796 };

   /**
    * Constructs a `BiNormalDist` object with default parameters @f$\mu_1
    * = \mu_2 = 0@f$, @f$\sigma_1 = \sigma_2 = 1@f$ and correlation
    * @f$\rho= @f$<tt> rho</tt>.
    */
   public BiNormalDist (double rho) {
      setParams (0.0, 1.0, 0.0, 1.0, rho);
   }

   /**
    * Constructs a `BiNormalDist` object with parameters @f$\mu_1@f$ =
    * `mu1`, @f$\mu_2@f$ = `mu2`, @f$\sigma_1@f$ = `sigma1`,
    * @f$\sigma_2@f$ = `sigma2` and @f$\rho@f$ = `rho`.
    */
   public BiNormalDist (double mu1, double sigma1,
                        double mu2, double sigma2, double rho) {
      setParams (mu1, sigma1, mu2, sigma2, rho);
   }


   public double density (double x, double y) {
      if (Math.abs(rho) == 1.)
         throw new IllegalArgumentException ("|rho| = 1");
      double X = (x - mu1)/sigma1;
      double Y = (y - mu2)/sigma2;
      double T = (X*X - 2.0*rho*X*Y + Y*Y) / (2.0*racRho*racRho);
      return Math.exp(-T) / detS;
   }

/**
 * Computes the *standard binormal* density function (
 * {@link REF_probdistmulti_BiNormalDist_eq_f1binormal
 * f1binormal} ) with @f$\mu_1 = \mu_2 = 0@f$ and @f$\sigma_1 = \sigma_2
 * = 1@f$.
 */
public static double density (double x, double y, double rho) {
       return density (0.0, 1.0, x, 0.0, 1.0, y, rho);
   }

   /**
    * Computes the *binormal* density function (
    * {@link REF_probdistmulti_BiNormalDist_eq_f1binormal
    * f1binormal} ) with parameters @f$\mu_1@f$ = `mu1`, @f$\mu_2@f$ =
    * `mu2`, @f$\sigma_1@f$ = `sigma1`, @f$\sigma_2@f$ = `sigma2` and
    * @f$\rho@f$ = `rho`.
    */
   public static double density (double mu1, double sigma1, double x, 
                                 double mu2, double sigma2, double y,
                                 double rho) {
      if (sigma1 <= 0.)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0.)
         throw new IllegalArgumentException ("sigma2 <= 0");
      if (Math.abs(rho) >= 1.)
         throw new IllegalArgumentException ("|rho| >= 1");
      double X = (x - mu1)/sigma1;
      double Y = (y - mu2)/sigma2;
      double fRho = (1.0 - rho)*(1.0 + rho);
      double T = (X*X - 2.0*rho*X*Y + Y*Y) / (2.0*fRho);
      return Math.exp (-T)/ (2.0*Math.PI*sigma1*sigma2*Math.sqrt(fRho));

   }


   protected static double Gauss (double z) {
       // Computes normal probabilities to within accuracy of 10^(-7)
       // Drezner Z., and G.O. Wesolowsky (1990) On the computation of the
       // bivariate normal integral.  J. Statist. Comput. Simul. 35:101-107.

       final double x = 1.0/(1.0 + 0.23164189 * Math.abs(z));
       double G = 0.53070271;
       for (int i = 0; i < 4; ++i)
          G = G*x + AGauss[i];
       G = G * x * Math.exp(-z*z/2.0);
       if (z > 0.0)
          G = 1.0 - G;
       return G;
   }


   protected static double specialCDF (double x, double y, double rho, double xbig) {
      // Compute the bivariate normal CDF for limiting cases and returns
      // its value. If non limiting case, returns -2 as flag.
      // xbig is practical infinity

      if (Math.abs (rho) > 1.0)
         throw new IllegalArgumentException ("|rho| > 1");
      if (x == 0.0 && y == 0.0)
         return 0.25 + Math.asin(rho)/(2.0*Math.PI);

      if (rho == 1.0) {
         if (y < x)
            x = y;
         return NormalDist.cdf01(x);
      }
      if (rho == -1.0) {
         if (y <= -x)
            return 0.0;
         else
            return NormalDist.cdf01(x) - NormalDist.cdf01(-y);
      }
      if (Math.abs (rho) < RHO_SMALL)
         // return NormalDist.cdf01(x) * NormalDist.cdf01(y);
         return Gauss(x) * Gauss(y); // return this value to keep sequence monotone around rho=0.

      if ((x <= -xbig) || (y <= -xbig))
         return 0.0;
      if (x >= xbig)
         return NormalDist.cdf01(y);
      if (y >= xbig)
         return NormalDist.cdf01(x);

      return -2.0;
   }

/**
 * Computes the standard *binormal* distribution (
 * {@link REF_probdistmulti_BiNormalDist_eq_cdf2binormal
 * cdf2binormal} ) using the fast Drezner-Wesolowsky method described in
 * @cite tDRE90a&thinsp;. The absolute error is expected to be smaller than
 * @f$2 \cdot10^{-7}@f$.
 */
public static double cdf (double x, double y, double rho) {
      double bvn = specialCDF (x, y, rho, 20.0);
      if (bvn >= 0.0)
         return bvn;
      bvn = 0.0;

      /* prob(x <= h1, y <= h2), where x and y are standard binormal, 
         with rho = corr(x,y),  error < 2e-7.
         Drezner Z., and G.O. Wesolowsky (1990) On the computation of the
         bivariate normal integral.  J. Statist. Comput. Simul. 35:101-107. */

      int i;
      double r1, r2, r3, rr, aa, ab, h3, h5, h6, h7;
      final double h1 = -x;
      double h2 = -y;
      final double h12 = (h1 * h1 + h2 * h2) / 2.;

      if (Math.abs (rho) >= 0.7) {
         r2 = (1. - rho) * (1. + rho);
         r3 = Math.sqrt (r2);
         if (rho < 0.)
            h2 = -h2;
         h3 = h1 * h2;
         if (h3 < 300.)
            h7 = Math.exp (-h3 / 2.);
         else
            h7 = 0.;
         if (r2 != 0.) {
            h6 = Math.abs (h1 - h2);
            h5 = h6 * h6 / 2.;
            h6 /= r3;
            aa = .5 - h3 / 8.;
            ab = 3. - 2. * aa * h5;
            bvn = .13298076 * h6 * ab * (1.0 - Gauss(h6))
               - Math.exp (-h5 / r2) * (ab + aa * r2) * 0.053051647;
//          if (h7 > 0.  && -h3 < 500.0)
            for (i = 0; i < 5; i++) {
               r1 = r3 * Z[i];
               rr = r1 * r1;
               r2 = Math.sqrt (1. - rr);
               bvn -= W[i] * Math.exp (-h5 / rr) * 
                   (Math.exp (-h3 / (1. + r2)) / r2 / h7 - 1. - aa * rr);
            }
         }

         if (rho > 0.)
            bvn = bvn * r3 * h7 + (1.0 - Gauss (Math.max (h1, h2)));
         else if (rho < 0.)
            bvn = (h1 < h2 ? Gauss (h2) - Gauss (h1) : 0.) - bvn * r3 * h7;

      } else {
         h3 = h1 * h2;
         for (i = 0; i < 5; i++) {
            r1 = rho * Z[i];
            r2 = 1. - r1 * r1;
            bvn += W[i] * Math.exp ((r1 * h3 - h12) / r2) / Math.sqrt (r2);
         }
         bvn = (1.0 - Gauss (h1)) * (1.0 - Gauss (h2)) + rho * bvn;
      }

      if (bvn <= 0.0)
         return 0.0;
      if (bvn <= 1.0)
         return bvn;
      return 1.0;
   }


   public double cdf (double x, double y) {
      return cdf ((x-mu1)/sigma1, (y-mu2)/sigma2, rho);
   }

/**
 * Computes the *binormal* distribution function (
 * {@link REF_probdistmulti_BiNormalDist_eq_cdf1binormal
 * cdf1binormal} ) with parameters @f$\mu_1@f$ = `mu1`, @f$\mu_2@f$ =
 * `mu2`, @f$\sigma_1@f$ = `sigma1`, @f$\sigma_2@f$ = `sigma2` and
 * @f$\rho@f$ = `rho`. Uses the fast Drezner-Wesolowsky method described in
 * @cite tDRE90a&thinsp;. The absolute error is expected to be smaller than
 * @f$2 \cdot10^{-7}@f$.
 */
public static double cdf (double mu1, double sigma1, double x, 
                             double mu2, double sigma2, double y,
                             double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      double X = (x - mu1)/sigma1;
      double Y = (y - mu2)/sigma2;
      return cdf (X, Y, rho);
   }

   /**
    * Computes the standard upper *binormal* distribution with @f$\mu_1 =
    * \mu_2 = 0@f$ and @f$\sigma_1 = \sigma_2 = 1@f$. Uses the fast
    * Drezner-Wesolowsky method described in @cite tDRE90a&thinsp;. The
    * absolute error is expected to be smaller than @f$2 \cdot10^{-7}@f$.
    */
   public static double barF (double x, double y, double rho) {
      return cdf (-x, -y, rho);
    }


   public double barF (double x, double y) {
      return barF ((x-mu1)/sigma1, (y-mu2)/sigma2, rho);
   }

/**
 * Computes the upper *binormal* distribution function (
 * {@link REF_probdistmulti_BiNormalDist_eq_cdf3binormal
 * cdf3binormal} ) with parameters @f$\mu_1@f$ = `mu1`, @f$\mu_2@f$ =
 * `mu2`, @f$\sigma_1@f$ = `sigma1`, @f$\sigma_2@f$ = `sigma2` and
 * @f$\rho@f$ = `rho`. Uses the fast Drezner-Wesolowsky method described in
 * @cite tDRE90a&thinsp;. The absolute error is expected to be smaller than
 * @f$2 \cdot10^{-7}@f$.
 */
public static double barF (double mu1, double sigma1, double x, 
                              double mu2, double sigma2, double y,
                              double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      double X = (x - mu1)/sigma1;
      double Y = (y - mu2)/sigma2;
      return barF (X, Y, rho);
   }


   public double[] getMean() {
      return getMean (mu1, mu2, sigma1, sigma2, rho);
   }

/**
 * Return the mean vector @f$E[X] = (\mu_1, \mu_2)@f$ of the binormal
 * distribution.
 */
public static double[] getMean(double mu1, double sigma1,
                                  double mu2, double sigma2, double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      if (Math.abs(rho) > 1.)
         throw new IllegalArgumentException ("|rho| > 1");

      double mean[] = new double[2];

      mean[0] = mu1;
      mean[1] = mu2;

      return mean;
   }


   public double[][] getCovariance() {
      return getCovariance (mu1, sigma1, mu2, sigma2, rho);
   }

/**
 * Return the covariance matrix of the binormal distribution.
 */
public static double[][] getCovariance (double mu1, double sigma1,
                                           double mu2, double sigma2,
                                           double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      if (Math.abs(rho) > 1.)
         throw new IllegalArgumentException ("|rho| > 1");

      double cov[][] = new double[2][2];

      cov[0][0] = sigma1 * sigma1;
      cov[0][1] = rho * sigma1 * sigma2;
      cov[1][0] = cov[0][1];
      cov[1][1] = sigma2 * sigma2;

      return cov;
   }


   public double[][] getCorrelation() {
      return getCovariance (mu1, sigma1, mu2, sigma2, rho);
   }

/**
 * Return the correlation matrix of the binormal distribution.
 */
public static double[][] getCorrelation (double mu1, double sigma1,
                                            double mu2, double sigma2,
                                            double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      if (Math.abs(rho) > 1.)
         throw new IllegalArgumentException ("|rho| > 1");

      double corr[][] = new double[2][2];

      corr[0][0] = 1.0;
      corr[0][1] = rho;
      corr[1][0] = rho;
      corr[1][1] = 1.0;

      return corr;
   }

   /**
    * Returns the parameter @f$\mu_1@f$.
    */
   public double getMu1() {
      return mu1;
   }

   /**
    * Returns the parameter @f$\mu_2@f$.
    */
   public double getMu2() {
      return mu2;
   }

   /**
    * Returns the parameter @f$\sigma_1@f$.
    */
   public double getSigma1() {
      return sigma1;
   }

   /**
    * Returns the parameter @f$\sigma_2@f$.
    */
   public double getSigma2() {
      return sigma2;
   }

   /**
    * Sets the parameters @f$\mu_1@f$ = `mu1`, @f$\mu_2@f$ = `mu2`,
    * @f$\sigma_1@f$ = `sigma1`, @f$\sigma_2@f$ = `sigma2` and
    * @f$\rho@f$ = `rho` of this object.
    */
   protected void setParams (double mu1, double sigma1, 
                             double mu2, double sigma2, double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      if (Math.abs(rho) > 1.)
         throw new IllegalArgumentException ("|rho| > 1");
      this.dimension = 2;
      this.mu1 = mu1;
      this.sigma1 = sigma1;
      this.mu2 = mu2;
      this.sigma2 = sigma2;
      this.rho = rho; 
      racRho = Math.sqrt((1.0 - rho)*(1.0 + rho));
      detS = 2.0*Math.PI*sigma1*sigma2*racRho;
   }
 
}