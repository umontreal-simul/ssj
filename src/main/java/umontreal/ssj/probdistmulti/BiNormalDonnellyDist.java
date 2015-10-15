/*
 * Class:        BiNormalDonnellyDist
 * Description:  bivariate normal distribution using Donnelly's code
 * Environment:  Java
 * Software:     SSJ 
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 */
package umontreal.ssj.probdistmulti;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.util.Num;

/**
 * Extends the class  @ref BiNormalDist for the *bivariate normal*
 * distribution @cite tJOH72a&thinsp; (page 84) using a translation of
 * Donnelly’s Fortran code in @cite tDON73a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_continuous2dim
 */
public class BiNormalDonnellyDist extends BiNormalDist {

private static final double TWOPI = 2.0*Math.PI;
private static final double SQRTPI = Math.sqrt(Math.PI);
private static final int KMAX = 6;

private static double BB[] = new double[KMAX + 1];
private static final double CB[] = {
    0.9999936, -0.9992989, 0.9872976,
   -0.9109973, 0.6829098, -0.3360210, 0.07612251 };


private static double BorthT (double h, double a)
{
   int k;
   final double w = a * h / Num.RAC2;
   final double w2 = w * w;
   double z = w * Math.exp(-w2);

   BB[0] = SQRTPI * (Gauss (Num.RAC2 * w) - 0.5);
   for (k = 1; k <= KMAX; ++k) {
      BB[k] = ((2 * k - 1) * BB[k - 1] - z) / 2.0;
      z *= w2;
   }

   final double h2 = h * h / 2.0;
   z = h / Num.RAC2;
   double sum = 0;
   for (k = 0; k <= KMAX; ++k) {
      sum += CB[k] * BB[k] / z;
      z *= h2;
   }
   return sum * Math.exp (-h2) / TWOPI;
}

   /**
    * Constructor with default parameters @f$\mu_1 = \mu_2 = 0@f$,
    * @f$\sigma_1 = \sigma_2 = 1@f$, correlation @f$\rho= @f$<tt>
    * rho</tt>, and @f$d = @f$ `ndig` digits of accuracy (the absolute
    * error is smaller than @f$10^{-d}@f$). Restriction: @f$d \le15@f$.
    */
   public BiNormalDonnellyDist (double rho, int ndig) {
       super (rho);
       if (ndig > 15)
          throw new IllegalArgumentException ("ndig > 15");
       decPrec = ndigit = ndig;
   }

   /**
    * Same as  {@link #BiNormalDonnellyDist(double,int)
    * BiNormalDonnellyDist(rho, 15)}.
    */
   public BiNormalDonnellyDist (double rho) {
       this (rho, 15);
   }

   /**
    * Constructor with parameters @f$\mu_1@f$ = `mu1`, @f$\mu_2@f$ =
    * `mu2`, @f$\sigma_1@f$ = `sigma1`, @f$\sigma_2@f$ = `sigma2`,
    * @f$\rho@f$ = `rho`, and @f$d = @f$ `ndig` digits of accuracy.
    * Restriction: @f$d \le15@f$.
    */
   public BiNormalDonnellyDist (double mu1, double sigma1, double mu2,
                                double sigma2, double rho, int ndig) {
      super (mu1, sigma1, mu2, sigma2, rho);
      if (ndig > 15)
         throw new IllegalArgumentException ("ndig > 15");
      decPrec = ndigit = ndig;
   }

   /**
    * Same as
    * {@link #BiNormalDonnellyDist(double,double,double,double,double,int)
    * BiNormalDonnellyDist(mu1, sigma1, mu2, sigma2, rho, 15)}.
    */
   public BiNormalDonnellyDist (double mu1, double sigma1, double mu2,
                                double sigma2, double rho) {
      this (mu1, sigma1, mu2, sigma2, rho, 15);
   }

   /**
    * The following methods use the parameter `ndig` for the number of
    * digits of absolute accuracy. If the same methods are called without
    * the `ndig` parameter, a default value of `ndig` = 15 will be used.
    *
    *  Computes the standard *binormal* distribution (
    * {@link REF_probdistmulti_BiNormalDist_eq_cdf2binormal
    * cdf2binormal} ) with the method described in @cite tDON73a&thinsp;,
    * where `ndig` is the number of decimal digits of accuracy provided
    * (<tt>ndig</tt> @f$\le15@f$). The code was translated from the
    * Fortran program written by T. G. Donnelly and copyrighted by the ACM
    * (see
    * [http://www.acm.org/pubs/copyright_policy/#Notice](http://www.acm.org/pubs/copyright_policy/#Notice)).
    * The absolute error is expected to be smaller than @f$10^{-d}@f$,
    * where @f$d=@f$ `ndig`.
    */
   public static double cdf (double x, double y, double rho, int ndig) {
   /* 
    * This is a translation of the FORTRAN code published by Thomas G. Donnelly
    * in the CACM, Vol. 16, Number 10, p. 638, (1973)
    */

      if (ndig > 15)
         throw new IllegalArgumentException ("ndig > 15");
      double b = specialCDF (x, y, rho, 13.0);
      if (b >= 0.0)
         return b;
      b = 0;

      final boolean SINGLE_FLAG = ndig <= 7 ? true : false;
      final double TWO_PI = 2.0 * Math.PI;
      final double r = rho;
      final double ah = -x;
      final double ak = -y;

      double a2, ap, cn, conex, ex, g2, gh, gk, gw = 0, h2, h4, rr, s1, s2,
         sgn, sn, sp, sqr, t, temp, w2, wh = 0, wk = 0;
      int is = -1;

      if (SINGLE_FLAG) {
         gh = Gauss (x) / 2.0;
         gk = Gauss (y) / 2.0;
      } else {
         gh = NormalDist.cdf01 (x) / 2.0;
         gk = NormalDist.cdf01 (y) / 2.0;
      }
      boolean flag = true;    // Easiest way to translate a Fortran goto

      rr = (1 - r) * (1 + r);
      if (rr < 0)
         throw new IllegalArgumentException ("|rho| > 1");
      sqr = Math.sqrt(rr);
      final double con = Math.PI * Num.TEN_NEG_POW[ndig];
      final double EPSILON = 0.5*Num.TEN_NEG_POW[ndig];

      if (ah != 0) {
         b = gh;
         if (ah * ak < 0)
            b -= 0.5;
         else if (ah * ak == 0) {
            flag = false;
         }
      } else if (ak == 0) {
         return Math.asin (r) / TWO_PI + 0.25;
      }

      if (flag)
         b += gk;
      if (ah != 0) {
         flag = false;
         wh = -ah;
         wk = (ak / ah - r) / sqr;
         gw = 2 * gh;
         is = -1;
      }

      do {
         if (flag) {
            wh = -ak;
            wk = (ah / ak - r) / sqr;
            gw = 2 * gk;
            is = 1;
         }
         flag = true;
         sgn = -1;
         t = 0;
         if (wk != 0) {
            if (Math.abs (wk) >= 1)
               if (Math.abs (wk) == 1) {
                  t = wk * gw * (1 - gw) / 2;
                  b += sgn * t;
                  if (is >= 0)        // Another Fortran goto
                     break;
                  else
                     continue;
               } else {
                  sgn = -sgn;
                  wh = wh * wk;
                  if (SINGLE_FLAG)
                     g2 = Gauss(wh);
                  else
                     g2 = NormalDist.cdf01(wh);
                  wk = 1 / wk;
                  if (wk < 0)
                     b = b + .5;
                  b = b - (gw + g2) / 2 + gw * g2;
               }
        /*****
            // Cette m'ethode de Borth est plus lente que simple Donnelly 
            if (ndig <= 7 && Math.abs (wh) > 1.6 && Math.abs (wk) > 0.3) {
               b += sgn * BorthT (wh, wk);
               if (is >= 0)
                  break;
               else
                  continue;
            }
        *****/
            h2 = wh * wh;
            a2 = wk * wk;
            h4 = h2 * .5;
            ex = 0;
            if (h4 < 300.0)
               ex = Math.exp (-h4);
            w2 = h4 * ex;
            ap = 1;
            s2 = ap - ex;
            sp = ap;
            s1 = 0;
            sn = s1;
            conex = Math.abs (con / wk);
            do {
               cn = ap * s2 / (sn + sp);
               s1 += cn;
               if (Math.abs (cn) <= conex)
                  break;
               sn = sp;
               sp += 1;
               s2 -= w2;
               w2 *= h4 / sp;
               ap = -ap * a2;
            } while (true);
            t = (Math.atan (wk) - wk * s1) / TWO_PI;
            b += sgn * t;
         }
         if (is >= 0)
            break;
      } while (ak != 0);

      if (b < EPSILON)
         b = 0;
      if (b > 1)
         b = 1;
      return b;
}

   /**
    * Computes the *binormal* distribution function (
    * {@link REF_probdistmulti_BiNormalDist_eq_cdf1binormal
    * cdf1binormal} ) with parameters @f$\mu_1@f$ = `mu1`, @f$\mu_2@f$ =
    * `mu2`, @f$\sigma_1@f$ = `sigma1`, @f$\sigma_2@f$ = `sigma2`,
    * correlation @f$\rho@f$ = `rho` and `ndig` decimal digits of
    * accuracy.
    */
   public static double cdf (double mu1, double sigma1, double x, 
                             double mu2, double sigma2, double y,
                             double rho, int ndig) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      double X = (x - mu1)/sigma1;
      double Y = (y - mu2)/sigma2;
      return cdf (X, Y, rho, ndig);
   }

   /**
    * Computes the upper *binormal* distribution function (
    * {@link REF_probdistmulti_BiNormalDist_eq_cdf3binormal
    * cdf3binormal} ) with parameters @f$\mu_1@f$ = `mu1`, @f$\mu_2@f$ =
    * `mu2`, @f$\sigma_1@f$ = `sigma1`, @f$\sigma_2@f$ = `sigma2`,
    * @f$\rho@f$ = `rho` and `ndig` decimal digits of accuracy.
    */
   public static double barF (double mu1, double sigma1, double x, 
                              double mu2, double sigma2, double y,
                              double rho, int ndig) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      double X = (x - mu1)/sigma1;
      double Y = (y - mu2)/sigma2;
      return barF (X, Y, rho, ndig);
   }

   /**
    * Computes the upper *standard binormal* distribution function (
    * {@link REF_probdistmulti_BiNormalDist_eq_cdf3binormal
    * cdf3binormal} ) with parameters @f$\rho@f$ = `rho` and `ndig`
    * decimal digits of accuracy.
    */
   public static double barF (double x, double y, double rho, int ndig)  {
      return cdf (-x, -y, rho, ndig);
   }
 

   public double cdf (double x, double y) {
      return cdf ((x-mu1)/sigma1, (y-mu2)/sigma2, rho, ndigit);
   }

   public static double cdf (double x, double y, double rho) {
       return cdf (x, y, rho, 15);
   }
 
   public static double cdf (double mu1, double sigma1, double x, 
                             double mu2, double sigma2, double y,
                             double rho) {
      return cdf (mu1, sigma1, x, mu2, sigma2, y, rho, 15);
   }

   public double barF (double x, double y) {
      return barF ((x-mu1)/sigma1, (y-mu2)/sigma2, rho, ndigit);
   }

   public static double barF (double mu1, double sigma1, double x, 
                              double mu2, double sigma2, double y,
                              double rho) {
      return barF (mu1, sigma1, x, mu2, sigma2, y, rho, 15);
   }

   public static double barF (double x, double y, double rho) {
      return barF (x, y, rho, 15);
   }
}