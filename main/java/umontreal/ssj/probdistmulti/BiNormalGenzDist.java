/*
 * Class:        BiNormalGenzDist
 * Description:  bivariate normal distribution using Genz's algorithm
 * Environment:  Java
 * Software:     SSJ
 * Organization: DIRO, Universite de Montreal
 * @author
 * @since
 */
package umontreal.ssj.probdistmulti;
import umontreal.ssj.probdist.NormalDist;

/**
 * Extends the class  @ref BiNormalDist for the *bivariate normal*
 * distribution @cite tJOH72a&thinsp; (page 84) using Genz’s algorithm as
 * described in @cite tGEN04a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_continuous2dim
 */
public class BiNormalGenzDist extends BiNormalDist {

   private static final double[][] W = {
   //       Gauss Legendre points and weights, n =  6
   { 0.1713244923791705, 0.3607615730481384, 0.4679139345726904},

   //       Gauss Legendre points and weights, n = 12
   { 0.4717533638651177e-1, 0.1069393259953183, 0.1600783285433464,
     0.2031674267230659, 0.2334925365383547, 0.2491470458134029 },

   //       Gauss Legendre points and weights, n = 20
   { 0.1761400713915212e-1, 0.4060142980038694e-1, 0.6267204833410906e-1,
     0.8327674157670475e-1, 0.1019301198172404, 0.1181945319615184,
     0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
     0.1527533871307259 }
   };

   private static final double[][] X = {
   //       Gauss Legendre points and weights, n =  6
   { 0.9324695142031522, 0.6612093864662647, 0.2386191860831970},

   //       Gauss Legendre points and weights, n = 12
   { 0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
     0.5873179542866171, 0.3678314989981802, 0.1252334085114692 },

   //       Gauss Legendre points and weights, n = 20
   { 0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
     0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
     0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
     0.7652652113349733e-1 }
   };

   /**
    * Constructs a `BiNormalGenzDist` object with default parameters
    * @f$\mu_1 = \mu_2 = 0@f$, @f$\sigma_1 = \sigma_2 = 1@f$ and
    * correlation @f$\rho= @f$<tt> rho</tt>.
    */
   public BiNormalGenzDist (double rho) {
       super (rho);
   }

   /**
    * Constructs a `BiNormalGenzDist` object with parameters @f$\mu_1@f$
    * = `mu1`, @f$\mu_2@f$ = `mu2`, @f$\sigma_1@f$ = `sigma1`,
    * @f$\sigma_2@f$ = `sigma2` and @f$\rho@f$ = `rho`.
    */
   public BiNormalGenzDist (double mu1, double sigma1,
                            double mu2, double sigma2, double rho) {
      super (mu1, sigma1, mu2, sigma2, rho);
   }

   /**
    * Computes the standard *binormal* distribution (
    * {@link REF_probdistmulti_BiNormalDist_eq_cdf2binormal
    * cdf2binormal} ) with the method described in @cite tGEN04a&thinsp;.
    * The code for the `cdf` was translated directly from the Matlab code
    * written by Alan Genz and available from his web page at
    * [http://www.math.wsu.edu/faculty/genz/homepage](http://www.math.wsu.edu/faculty/genz/homepage)
    * (the code is copyrighted by Alan Genz and is included in this
    * package with the kind permission of the author). The absolute error
    * is expected to be smaller than @f$0.5 \cdot10^{-15}@f$.
    */
   public static double cdf (double x, double y, double rho) {
      double bvn = specialCDF (x, y, rho, 40.0);
      if (bvn >= 0.0)
         return bvn;

/*
//   Copyright (C) 2005, Alan Genz,  All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided the following conditions are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The contributor name(s) may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
//   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
//   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
//   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//   function p = bvnl( dh, dk, r )
//
//  A function for computing bivariate normal probabilities.
//  bvnl calculates the probability that x < dh and y < dk.
//    parameters
//      dh 1st upper integration limit
//      dk 2nd upper integration limit
//      r   correlation coefficient
//
//   Author
//       Alan Genz
//       Department of Mathematics
//       Washington State University
//       Pullman, Wa 99164-3113
//       Email : alangenz@wsu.edu
//   This function is based on the method described by
//        Drezner, Z and G.O. Wesolowsky, (1989),
//        On the computation of the bivariate normal inegral,
//        Journal of Statist. Comput. Simul. 35, pp. 101-107,
//    with major modifications for double precision, for |r| close to 1,
//    and for matlab by Alan Genz - last modifications 7/98.
//
//      p = bvnu( -dh, -dk, r );
//      return
//
//   end bvnl
//
//      function p = bvnu( dh, dk, r )
//
//  A function for computing bivariate normal probabilities.
//  bvnu calculates the probability that x > dh and y > dk.
//    parameters
//      dh 1st lower integration limit
//      dk 2nd lower integration limit
//      r   correlation coefficient
//
//   Author
//       Alan Genz
//       Department of Mathematics
//       Washington State University
//       Pullman, Wa 99164-3113
//       Email : alangenz@wsu.edu
//
//    This function is based on the method described by
//        Drezner, Z and G.O. Wesolowsky, (1989),
//        On the computation of the bivariate normal inegral,
//        Journal of Statist. Comput. Simul. 35, pp. 101-107,
//    with major modifications for double precision, for |r| close to 1,
//    and for matlab by Alan Genz - last modifications 7/98.
//        Note: to compute the probability that x < dh and y < dk, use
//              bvnu( -dh, -dk, r ).
//
*/

      final double TWOPI = 2.0 * Math.PI;
      final double sqrt2pi = 2.50662827463100050241; // sqrt(TWOPI)
      double h, k, hk, hs, asr, sn, as, a, b, c, d, sp, rs, ep, bs, xs;
      int i, lg, ng, is;

      if (Math.abs (rho) < 0.3) {
         ng = 0;
         lg = 3;

      } else if (Math.abs (rho) < 0.75) {
         ng = 1;
         lg = 6;

      } else {
         ng = 2;
         lg = 10;
      }

      h = -x;
      k = -y;
      hk = h * k;
      bvn = 0;
      if (Math.abs (rho) < 0.925) {
         hs = (h * h + k * k) / 2.0;
         asr = Math.asin (rho);
         for (i = 0; i < lg; ++i) {
            sn = Math.sin (asr * (1.0 - X[ng][i]) / 2.0);
            bvn += W[ng][i] * Math.exp ((sn * hk - hs) / (1.0 - sn * sn));
            sn = Math.sin (asr * (1.0 + X[ng][i]) / 2.0);
            bvn += W[ng][i] * Math.exp ((sn * hk - hs) / (1.0 - sn * sn));
         }
         bvn =  bvn * asr /(4.0*Math.PI) +
                NormalDist.cdf01 (-h) * NormalDist.cdf01 (-k);

      } else {
         if (rho < 0.0) {
            k = -k;
            hk = -hk;
         }
         if (Math.abs (rho) < 1.0) {
            as = (1.0 - rho) * (1.0 + rho);
            a = Math.sqrt (as);
            bs = (h - k) * (h - k);
            c = (4.0 - hk) / 8.0;
            d = (12.0 - hk) / 16.0;
            asr = -(bs / as + hk) / 2.0;
            if (asr > -100.0)
               bvn = a * Math.exp (asr) * (1.0 - c * (bs - as) * (1.0 -
                     d * bs / 5.0) / 3.0 + c * d * as * as / 5.0);

            if (-hk < 100.0) {
               b = Math.sqrt (bs);
               sp = sqrt2pi * NormalDist.cdf01 (-b / a);
               bvn = bvn - Math.exp (-hk / 2.0) * sp * b * (1.0 - c * bs * (1.0 -
                     d * bs / 5.0) / 3.0);
            }
            a = a / 2.0;
            for (i = 0; i < lg; ++i) {
               for (is = -1; is <= 1; is += 2) {
                  xs = (a * (is * X[ng][i] + 1.0));
                  xs = xs * xs;
                  rs = Math.sqrt (1.0 - xs);
                  asr = -(bs / xs + hk) / 2.0;
                  if (asr > -100.0) {
                     sp = (1.0 + c * xs * (1.0 + d * xs));
                     ep = Math.exp (-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs;
                     bvn += a * W[ng][i] * Math.exp (asr) * (ep - sp);
                  }
               }
            }
            bvn = -bvn / TWOPI;
         }
         if (rho > 0.0) {
            if (k > h)
               h = k;
            bvn += NormalDist.cdf01 (-h);
         }
         if (rho < 0.0) {
            xs = NormalDist.cdf01(-h) - NormalDist.cdf01(-k);
            if (xs < 0.0)
               xs = 0.0;
            bvn = -bvn + xs;
         }
      }
   if (bvn <= 0.0)
      return 0.0;
   if (bvn >= 1.0)
      return 1.0;
   return bvn;

   }


   public static double cdf (double mu1, double sigma1, double x,
                             double mu2, double sigma2, double y,
                             double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      double Z = (x - mu1)/sigma1;
      double T = (y - mu2)/sigma2;
      return cdf (Z, T, rho);
   }

   public double cdf (double x, double y) {
      return cdf ((x-mu1)/sigma1, (y-mu2)/sigma2, rho);
   }

   public double barF (double x, double y) {
      return barF ((x-mu1)/sigma1, (y-mu2)/sigma2, rho);
   }

   public static double barF (double mu1, double sigma1, double x,
                              double mu2, double sigma2, double y,
                              double rho) {
      if (sigma1 <= 0)
         throw new IllegalArgumentException ("sigma1 <= 0");
      if (sigma2 <= 0)
         throw new IllegalArgumentException ("sigma2 <= 0");
      double Z = (x - mu1)/sigma1;
      double T = (y - mu2)/sigma2;
      return barF (Z, T, rho);
   }

   public static double barF (double x, double y, double rho) {
      return cdf (-x, -y, rho);
   }
}