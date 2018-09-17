/*
 * Class:        BetaDist
 * Description:  beta distribution
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
package  umontreal.ssj.probdist;
import umontreal.ssj.util.*;
import optimization.*;

/**
 * Extends the class  @ref ContinuousDistribution for the *beta* distribution
 * @cite tJOH95b&thinsp; (page 210) with shape parameters @f$\alpha> 0@f$
 * and @f$\beta> 0@f$, over the interval @f$[a,b]@f$, where @f$a < b@f$.
 * This distribution has density
 * @f[
 *   f(x) = \frac{ (x-a)^{\alpha- 1}(b - x)^{\beta- 1}}{\mathcal{B} (\alpha, \beta)(b - a)^{\alpha+ \beta- 1}}, \qquad\mbox{for } a\le x\le b, \mbox{ and }0\mbox{ elsewhere},
 * @f]
 * and distribution function
 * @anchor REF_probdist_BetaDist_eq_Fbeta
 * @f[
 *   F(x) = I_{\alpha,\beta}(x) = \int_a^x \frac{(\xi- a)^{\alpha-1} (b - \xi)^{\beta-1}}{\mathcal{B} (\alpha, \beta)(b - a)^{\alpha+ \beta- 1}} d\xi, \qquad\mbox{for } a\le x\le b, \tag{Fbeta}
 * @f]
 * where @f$\mathcal{B}(\alpha,\beta)@f$ is the *beta* function defined by
 * @anchor REF_probdist_BetaDist_eq_betadef
 * @f[
 *   \mathcal{B} (\alpha,\beta) = \frac{\Gamma(\alpha) \Gamma(\beta)}{ \Gamma(\alpha+\beta)},\tag{betadef}
 * @f]
 * and @f$\Gamma(x)@f$ is the gamma function defined in (
 * {@link REF_probdist_GammaDist_eq_Gamma Gamma} ).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class BetaDist extends ContinuousDistribution {
   protected double alpha;         // First parameter
   protected double beta;          // Second parameter
   protected double a, b;          // Interval x in [a, b]
   protected double bminusa;
   protected double logFactor;
   protected double Beta;          // Function Beta(alpha, beta)
   protected double logBeta;       // Ln(Beta(alpha, beta))

   private static class Optim implements Lmder_fcn
   {
      private double a;
      private double b;

      public Optim (double a, double b)
      {
         this.a = a;
         this.b = b;
      }

      public void fcn (int m, int n, double[] x, double[] fvec, double[][] fjac, int iflag[])
      {
         if (x[1] <= 0.0 || x[2] <= 0.0) {
             final double BIG = 1.0e100;
             fvec[1] = BIG;
             fvec[2] = BIG;
             fjac[1][1] = BIG;
             fjac[1][2] = 0.0;
             fjac[2][1] = 0.0;
             fjac[2][2] = BIG;
             return;
         }

         double trig;
         if (iflag[1] == 1)
         {
            trig = Num.digamma (x[1] + x[2]);
            fvec[1] = Num.digamma(x[1]) - trig - a;
            fvec[2] = Num.digamma(x[2]) - trig - b;
         }
         else if (iflag[1] == 2)
         {
            trig = Num.trigamma (x[1] + x[2]);

            fjac[1][1] = Num.trigamma (x[1]) - trig;
            fjac[1][2] = - trig;
            fjac[2][1] = - trig;
            fjac[2][2] = Num.trigamma (x[2]) - trig;
         }
      }
   }

   /**
    * Constructs a `BetaDist` object with parameters @f$\alpha=@f$
    * `alpha`, @f$\beta=@f$ `beta` and default domain @f$[0,1]@f$.
    */
   public BetaDist (double alpha, double beta) {
      setParams (alpha, beta, 0.0, 1.0);
   }

   /**
    * Constructs a `BetaDist` object with parameters @f$\alpha=@f$
    * `alpha`, @f$\beta=@f$ `beta` and domain
    * @f$[@f$<tt>a</tt>@f$,@f$&nbsp;<tt>b</tt>@f$]@f$.
    */
   public BetaDist (double alpha, double beta, double a, double b) {
      setParams (alpha, beta, a, b);
   }
   @Deprecated
   public BetaDist (double alpha, double beta, int d) {
      setParams (alpha, beta, 0.0, 1.0, d);
   }

   @Deprecated
   public BetaDist (double alpha, double beta, double a, double b, int d) {
      setParams (alpha, beta, a, b, d);
   }


   @Override
   public double density (double x) {
      if (x <= a || x >= b)
         return 0;
      double temp = (alpha - 1) * Math.log(x - a) + (beta - 1) * Math.log(b - x);
//      return factor*Math.pow (x - a, alpha - 1)*Math.pow (b - x, beta - 1);
      return Math.exp(logFactor + temp);
   }

   @Override
   public double cdf (double x) {
      return cdf (alpha, beta, (x - a)/bminusa);
   }

   @Override
   public double barF (double x) {
      return barF (alpha, beta, (x - a)/bminusa);
   }

   @Override
   public double inverseF (double u) {
      return a + (b - a)*inverseF (alpha, beta, u);
   }

   @Override
   public double getMean() {
      return BetaDist.getMean (alpha, beta, a, b);
   }

   @Override
   public double getVariance() {
      return BetaDist.getVariance (alpha, beta, a, b);
   }

   @Override
   public double getStandardDeviation() {
      return BetaDist.getStandardDeviation (alpha, beta, a, b);
   }

/**
 * Same as  {@link #density(double,double,double,double,double)
 * density(alpha, beta, 0, 1, x)}.
 */
public static double density (double alpha, double beta, double x) {
      return density (alpha, beta, 0.0, 1.0, x);
   }

   /**
    * Computes the density function of the *beta* distribution.
    */
   public static double density (double alpha, double beta,
                                 double a, double b, double x) {
      if (a >= b)
         throw new IllegalArgumentException ("a >= b");
      if (x <= a || x >= b)
         return 0;

      double z = -Num.lnBeta (alpha, beta) - (alpha + beta - 1)* Math.log(b-a) +
		      (alpha-1)*Math.log(x-a) + (beta-1)*Math.log(b-x);
      return Math.exp(z);
   }


   static double beta_g (double x)
   /*
    * Used in the normal approximation of beta. This is the function
    * (1 - x^2 + 2x*ln(x)) / (1 - x)^2.
    */
   {
      if (x > 1.0)
         return -beta_g(1.0/x);
      if (x < 1.0e-200)
         return 1.0;

      final double y = 1.0 - x;
      if (x < 0.9)
         return (1.0 - x*x + 2.0*x*Math.log(x)) / (y*y);
      if (x == 1.0)
         return 0.0;

      // For x near 1, use a series expansion to avoid loss of precision
      final double EPS = 1.0e-12;
      double term;
      double ypow = 1.0;
      double sum = 0.0;
      int j = 2;
      do {
         ypow *= y;
         term = ypow / (j * (j + 1));
         sum += term;
         j++;
      } while (Math.abs (term / sum) > EPS);

      return 2.0 * sum;
   }


   private static double bolshev (double alpha, double beta, int d, double x) {
      // Bol'shev approximation for large max (alpha, beta)
      // and small min (alpha, beta)
     /* if (x > 0.5)
         return barF (beta, alpha, 1.0 - x); */
      boolean flag = false;
      double u, temp, yd, gam;

      if (alpha < beta) {
         u = alpha;
         alpha = beta;
         beta = u;
         flag = false;
      } else
         flag = true;

      u = alpha + 0.5 * beta - 0.5;
      if (!flag)
         temp = x / (2.0 - x);
      else
         temp = (1.0 - x) / (1.0 + x);
      yd = 2.0 * u * temp;
      gam = (Math.exp (beta * Math.log (yd) - yd -
             Num.lnGamma (beta)) * (2.0 * yd * yd - (beta - 1.0) * yd -
             (beta * beta - 1.0))) / (24.0 * u * u);
      if (flag) {
         yd = GammaDist.barF (beta, d, yd);
         return Math.max(0, yd - gam);
      } else {
         yd = GammaDist.cdf (beta, d, yd);
         return yd + gam;
      }
   }


   private static double peizer (double alpha, double beta, double x) {
      // Normal approximation of Peizer and Pratt.   Reference: @cite tPEI68a
      double temp, h1, h3, y;
      h1 = alpha + beta - 1.0;
      y = 1.0 - x;
      if (x > 1.0e-15)
         temp = 1.0 + y*beta_g ((alpha - 0.5)/(h1*x));
      else
         temp = GammaDist.mybelog ((alpha - 0.5) / (h1 * x));

      h3 = Math.sqrt ((temp + x*beta_g ((beta - 0.5)/(h1*y)))
         /((h1 + 1.0/6.0)*x*y))
         *((h1 + 1.0/3.0 + 0.02*(1.0/alpha + 1.0/beta + 1.0/(alpha + beta)))
         *x - alpha + 1.0/3.0 - 0.02/alpha - 0.01/(alpha + beta));

      return NormalDist.cdf01 (h3);
   }


   private static double donato (double alpha, double beta, double x) {
      // Cuyt p. 387, 18.5.20b
      // distribution Beta avec fractions continues
      // Il faut choisir MMAX >= sqrt(max(alpha, beta))

      double mid = (alpha + 1.0) / (alpha + beta + 2.0);
      if (x > mid)
         return 1.0 - donato (beta, alpha, 1.0 - x);

      final int MMAX = 100;    // pour ALPHABETAMAX = 10000
      double[] Ta = new double[1 + MMAX];
      double[] Tb = new double[1 + MMAX];
      int M1 = MMAX;
      double tem, tem1;
      int m;

      if ((beta <= MMAX) && (beta % 1.0 < 1.0e-100)) {
         // if beta is integer, Ta[i0] = 0; it is useless to evaluate
         // the other Ta[i] for i > i0
         M1 = (int) beta;
      }

      Ta[1] = 1;
      for (m = 1; m < M1; m++) {
         tem = alpha + 2 * m - 1;
         Ta[m + 1] = (alpha + m - 1) * (alpha + beta + m - 1) *
                     (beta - m) * m * x * x / (tem * tem);
      }

      // term m = 0 in the next loop; avoid tem1 = 0/0 for alpha = 1
      tem = alpha * (alpha + beta) / (alpha + 1);
      Tb[1] = alpha - tem * x;

      for (m = 1; m < M1; m++) {
         tem = (alpha + m) * (alpha + beta + m) / (alpha + 2 * m + 1);
         tem1 = m * (beta - m) / (alpha + 2 * m - 1);
         Tb[m + 1] = alpha + 2 * m + (tem1 - tem) * x;
      }

      while (0. == Tb[M1] && M1 > 1) {
         --M1;
      }

      // evaluate continuous fraction
      double con = 0;
      for (m = M1; m > 0; m--) {
         con = Ta[m] / (Tb[m] + con);
      }

      tem = Num.lnBeta (alpha, beta) - alpha * Math.log (x) - beta * Math.log1p (-x);
      return con * Math.exp (-tem);
   }

   @Deprecated
   public static double cdf (double alpha, double beta, int d, double x) {
      return cdf (alpha, beta, x);
   }


   @Deprecated
   public static double cdf (double alpha, double beta,
                             double a, double b, int d, double x) {
      return cdf (alpha, beta, d, (x - a)/(b - a));
   }


   @Deprecated
   public static double barF (double alpha, double beta, int d, double x) {
      return 1.0 - cdf (alpha, beta, d, x);
   }


   @Deprecated
   public static double barF (double alpha, double beta,
                              double a, double b, int d, double x) {
      if (a >= b)
         throw new IllegalArgumentException ("a >= b");
      return 1.0 - cdf (alpha, beta, d, (x - a)/(b - a));
   }

/**
 * Same as  {@link #cdf(double,double,double,double,double) cdf(alpha, beta,
 * 0, 1, x)}.
 */
public static double cdf (double alpha, double beta, double x) {
      if (alpha <= 0.0)
        throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
        throw new IllegalArgumentException ("beta <= 0");

      if (x <= 0.0)
         return 0.0;
      if (x >= 1.0)
         return 1.0;
      if (1.0 == beta)
         return Math.pow(x, alpha);

      final double ALPHABETAMAX = 10000.0;
      final double ALPHABETALIM = 30.0;

      if (Math.max (alpha, beta) <= ALPHABETAMAX) {
         return donato (alpha, beta, x);
      }

      if ((alpha > ALPHABETAMAX && beta < ALPHABETALIM) ||
          (beta > ALPHABETAMAX && alpha < ALPHABETALIM)) {
         return bolshev (alpha, beta, 12, x);
      }

      return peizer (alpha, beta, x);
   }

   /**
    * Computes the distribution function.  If @f$\max(\alpha, \beta)
    * \le10^4@f$, uses a continuous fraction in @f$\alpha@f$ and
    * @f$\beta@f$ given in @cite tDID92a, @cite tCUY08a&thinsp;.
    * Otherwise, if @f$\min(\alpha, \beta) \le30@f$, uses an
    * approximation due to Bol’shev @cite tMAR78a&thinsp;, else uses a
    * normal approximation @cite tPEI68a&thinsp;.
    */
   public static double cdf (double alpha, double beta,
                             double a, double b, double x) {
      return cdf (alpha, beta, (x - a)/(b - a));
   }

   /**
    * Same as  {@link #barF(double,double,double,double,double)
    * barF(alpha, beta, 0, 1, x)}.
    */
   public static double barF (double alpha, double beta, double x) {
      return cdf (beta, alpha, 1.0 - x);
   }

   /**
    * Computes the complementary distribution function.
    */
   public static double barF (double alpha, double beta,
                              double a, double b, double x) {
      if (a >= b)
         throw new IllegalArgumentException ("a >= b");
      return cdf (beta, alpha, (b - x)/(b - a));
   }


   @Deprecated
   public static double inverseF (double alpha, double beta, int d, double u) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (d <= 0)
        throw new IllegalArgumentException ("d <= 0");
      if (u > 1.0 || u < 0.0)
         throw new IllegalArgumentException ("u not in [0,1]");
      if (u <= 0)
         return 0;
      if (u >= 1)
         return 1;

      /*
       * Code taken from
       * Cephes Math Library Release 2.8:  June, 2000
       * Copyright 1984, 1996, 2000 by Stephen L. Moshier
       */
      final double MACHEP = 1.11022302462515654042E-16;
      final double MAXLOG = 7.09782712893383996843E2;
      final double MINLOG = -7.08396418532264106224E2;
    //  final double MAXNUM = 1.7976931348623158E308;

      boolean ihalve = false;
      boolean newt = false;

      double p = 0, q = 0, y0 = 0, z = 0, y = 0, x = 0, x0, x1, lgm = 0,
             yp = 0, di = 0, dithresh = 0, yl, yh, xt = 0;
      int i, dir;
      boolean rflg, nflg;
      x0 = 0.0;
      yl = 0.0;
      x1 = 1.0;
      yh = 1.0;
      nflg = false;
      rflg = false;
      if (alpha <= 1.0 || beta <= 1.0) {
         dithresh = 1.0e-6;
         rflg = false;
         p = alpha;
         q = beta;
         y0 = u;
         x = p/(p+q);
         y = cdf (p, q, x);
         ihalve = true;
      }
      else
         dithresh = 1.0e-4;

mainloop:
      while (true) {
         if (ihalve) {
            ihalve = false;
            dir = 0;
            di = 0.5;
            for (i = 0; i<100; i++) {
               if (i != 0) {
                  x = x0  +  di*(x1 - x0);
                  if (x == 1.0)
                     x = 1.0 - MACHEP;
                  if (x == 0.0) {
                     di = 0.5;
                     x = x0  +  di*(x1 - x0);
                     if (x == 0.0) {
                        // System.err.println ("BetaDist.inverseF: underflow");
                        return 0;
                     }
                  }
                  y = cdf (p, q, x);
                  yp = (x1 - x0)/(x1 + x0);
                  if (Math.abs (yp) < dithresh) {
                     newt = true;
                     continue mainloop;
                  }
                  yp = (y-y0)/y0;
                  if (Math.abs (yp) < dithresh) {
                     newt = true;
                     continue mainloop;
                  }
               }
               if (y < y0) {
                  x0 = x;
                  yl = y;
                  if (dir < 0) {
                     dir = 0;
                     di = 0.5;
                  }
                  else if (dir > 3)
                     di = 1.0 - (1.0 - di)*(1.0 - di);
                  else if (dir > 1)
                     di = 0.5*di + 0.5;
                  else
                     di = (y0 - y)/(yh - yl);
                  dir += 1;
                  if (x0 > 0.75) {
                  // if (0 == y)
                  //    y = EPS;
                     if (rflg) {
                        rflg = false;
                        p = alpha;
                        q = beta;
                        y0 = u;
                     }
                     else {
                        rflg = true;
                        p = beta;
                        q = alpha;
                        y0 = 1.0 - u;
                     }
                     x = 1.0 - x;
                     y = cdf (p, q, x);
                     x0 = 0.0;
                     yl = 0.0;
                     x1 = 1.0;
                     yh = 1.0;
                     ihalve = true;
                     continue mainloop;
                  }
               }
               else {
                  x1 = x;
                  if (rflg && x1 < MACHEP) {
                     x = 0.0;
                     break mainloop;
                  }
                  yh = y;
                  if (dir > 0) {
                     dir = 0;
                     di = 0.5;
                  }
                  else if (dir < -3)
                     di = di*di;
                  else if (dir < -1)
                     di = 0.5*di;
                  else
                     di = (y - y0)/(yh - yl);
                  dir -= 1;
               }
            }
            // PLOSS error
            if (x0 >= 1.0) {
               x = 1.0 - MACHEP;
               break mainloop;
            }
            if (x <= 0.0) {
            // System.err.println ("BetaDist.inverseF: underflow");
               return 0 ;
            }
            newt = true;
         }
         if (newt) {
            newt = false;
            if (nflg)
               break mainloop;
            nflg = true;
            lgm = Num.lnGamma (p+q) - Num.lnGamma (p) - Num.lnGamma (q);

            for (i=0; i<8; i++) {
               /* Compute the function at this point. */
               if (i != 0)
                  y = cdf (p, q, x);
               if (y < yl) {
                  x = x0;
                  y = yl;
               }
               else if (y > yh) {
                  x = x1;
                  y = yh;
               }
               else if (y < y0) {
                  x0 = x;
                  yl = y;
               }
               else {
                  x1 = x;
                  yh = y;
               }
               if (x >= 1.0 || x <= 0.0)
                  break;
               /* Compute the derivative of the function at this point. */
               z = (p - 1.0)*Math.log (x) + (q - 1.0)*Math.log1p (-x) + lgm;
               if (z < MINLOG)
                  break mainloop;
               if (z > MAXLOG)
                  break;
               z = Math.exp (z);
               /* Compute the step to the next approximation of x. */
               z = (y - y0)/z;
               xt = x - z;
               if (xt <= x0) {
                  y = (x - x0) / (x1 - x0);
                  xt = x0 + 0.5*y*(x - x0);
                  if (xt <= 0.0)
                     break;
               }
               if (xt >= x1) {
                  y = (x1 - x) / (x1 - x0);
                  xt = x1 - 0.5*y*(x1 - x);
                  if (xt >= 1.0)
                     break;
               }
               x = xt;
               if (Math.abs (z/x) < 128.0*MACHEP)
                  break mainloop;
            }
            /* Did not converge.  */
            dithresh = 256.0*MACHEP;
            ihalve = true;
            continue mainloop;
         }

         yp = -NormalDist.inverseF01 (u);

         if (u > 0.5) {
            rflg = true;
            p = beta;
            q = alpha;
            y0 = 1.0 - u;
            yp = -yp;
         }
         else {
            rflg = false;
            p = alpha;
            q = beta;
            y0 = u;
         }

         lgm = (yp*yp - 3.0)/6.0;
         x = 2.0/(1.0/(2.0*p-1.0)  +  1.0/(2.0*q-1.0));
         z = yp*Math.sqrt (x + lgm)/x
           - (1.0/(2.0*q-1.0) - 1.0/(2.0*p-1.0) )
           * (lgm + 5.0/6.0 - 2.0/(3.0*x));
         z = 2.0*z;
         if (z < MINLOG) {
            x = 1.0;
            // System.err.println ("BetaDist.inverseF: underflow");
            return 0;
         }
         x = p/( p + q*Math.exp (z));
         y = cdf (p, q, x);
         yp = (y - y0)/y0;
         if (Math.abs (yp) < 0.2) {
            newt = true;
            continue mainloop;
         }
         ihalve = true;
      }

      // Done
      if (rflg) {
         if (x <= MACHEP)
            x = 1.0 - MACHEP;
         else
            x = 1.0 - x;
      }
      return x;
   }

/**
 * Same as  {@link #inverseF(double,double,double,double,double)
 * inverseF(alpha, beta, 0, 1, u)}.
 */
public static double inverseF (double alpha, double beta, double u) {
      return inverseF (alpha, beta, 12, u);
   }


   @Deprecated
   public static double inverseF (double alpha, double beta,
                                  double a, double b, int d, double u) {
      if (a >= b)
        throw new IllegalArgumentException ("a >= b");
      return a + (b - a)*inverseF (alpha, beta, d, u);
   }

/**
 * Returns the inverse beta distribution function using the algorithm
 * implemented in @cite iMOS00a&thinsp;. The method performs interval halving
 * or Newton iterations to compute the inverse.
 */
public static double inverseF (double alpha, double beta,
                                  double a, double b, double u) {
      if (a >= b)
        throw new IllegalArgumentException ("a >= b");
      return a + (b - a)*inverseF (alpha, beta, u);
   }

   /**
    * Estimates the parameters @f$(\alpha,\beta)@f$ of the beta
    * distribution over the interval @f$[0,1]@f$ using the maximum
    * likelihood method, from the @f$n@f$ observations @f$x[i]@f$, @f$i =
    * 0, 1,…, n-1@f$. The estimates are returned in a two-element array,
    * in regular order: [@f$\alpha@f$, @f$\beta@f$].  The maximum
    * likelihood estimators are the values @f$(\hat{\alpha},
    * \hat{\beta})@f$ that satisfy the equations:
    * @f{align*}{
    *    \psi(\alpha) - \psi(\alpha+ \beta) 
    *    & 
    *    = 
    *    \frac{1}{n} \sum_{i=1}^n \ln(x_i)
    *    \\ 
    *   \psi(\beta) - \psi(\alpha+ \beta) 
    *    & 
    *    = 
    *    \frac{1}{n} \sum_{i=1}^n \ln(1 - x_i)
    * @f}
    * where @f$\bar{x}_n@f$ is the average of @f$x[0],…,x[n-1]@f$, and
    * @f$\psi@f$ is the logarithmic derivative of the Gamma function
    * @f$\psi(x) = \Gamma’(x) / \Gamma(x)@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    *  @return returns the parameters [@f$\hat{\alpha}@f$,
    * @f$\hat{\beta}@f$]
    */
   public static double[] getMLE (double[] x, int n) {
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");

      double sum = 0.0;
      double a = 0.0;
      double b = 0.0;
      for (int i = 0; i < n; i++)
      {
         sum += x[i];
         if (x[i] > 0.0)
            a += Math.log (x[i]);
         else
            a -= 709.0;
         if (x[i] < 1.0)
            b += Math.log1p (-x[i]);
         else
            b -= 709.0;
      }
      double mean = sum / n;

      sum = 0.0;
      for (int i = 0; i < n; i++)
         sum += (x[i] - mean) * (x[i] - mean);
      double var = sum / (n - 1);

      Optim system = new Optim (a, b);

      double[] param = new double[3];
      param[1] = mean * ((mean * (1.0 - mean) / var) - 1.0);
      param[2] = (1.0 - mean) * ((mean * (1.0 - mean) / var) - 1.0);
      double[] fvec = new double [3];
      double[][] fjac = new double[3][3];
      int[] info = new int[2];
      int[] ipvt = new int[3];

      Minpack_f77.lmder1_f77 (system, 2, 2, param, fvec, fjac, 1e-5, info, ipvt);

      double parameters[] = new double[2];
      parameters[0] = param[1];
      parameters[1] = param[2];

      return parameters;
   }

   /**
    * Creates a new instance of a beta distribution with parameters
    * @f$\alpha@f$ and @f$\beta@f$ over the interval @f$[0,1]@f$
    * estimated using the maximum likelihood method based on the @f$n@f$
    * observations @f$x[i]@f$, @f$i = 0, 1, …, n-1@f$.
    *  @param x            the list of observations to use to evaluate
    *                      parameters
    *  @param n            the number of observations to use to evaluate
    *                      parameters
    */
   public static BetaDist getInstanceFromMLE (double[] x, int n) {
      double parameters[] = getMLE (x, n);
      return new BetaDist (parameters[0], parameters[1]);
   }

   /**
    * Computes and returns the mean @f$E[X] = \alpha/ (\alpha+
    * \beta)@f$ of the beta distribution with parameters @f$\alpha@f$
    * and @f$\beta@f$, over the interval @f$[0, 1]@f$.
    *  @return the mean of the Beta distribution
    */
   public static double getMean (double alpha, double beta) {
      return getMean (alpha, beta, 0.0, 1.0);
   }

   /**
    * Computes and returns the mean @f$E[X] = (b\alpha+ a\beta)/
    * (\alpha+ \beta)@f$ of the beta distribution with parameters
    * @f$\alpha@f$ and @f$\beta@f$ over the interval @f$[a, b]@f$.
    *  @return the mean of the Beta distribution
    */
   public static double getMean (double alpha, double beta, double a,
                                 double b) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      return (alpha*b + beta*a) / (alpha + beta);
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] =
    * \frac{\alpha\beta}{(\alpha+ \beta)^2 (\alpha+ \beta+ 1)}@f$ of
    * the beta distribution with parameters @f$\alpha@f$ and
    * @f$\beta@f$, over the interval @f$[0, 1]@f$.
    *  @return the variance of the beta distribution @f$\mbox{Var}[X] =
    * \alpha\beta/ [(\alpha+ \beta)^2 (\alpha+ \beta+ 1)]@f$.
    */
   public static double getVariance (double alpha, double beta) {
      return getVariance (alpha, beta, 0.0, 1.0);
   }

   /**
    * Computes and returns the variance @f$\mbox{Var}[X] =
    * \frac{\alpha\beta(b-a)^2}{(\alpha+ \beta)^2 (\alpha+ \beta+
    * 1)}@f$ of the beta distribution with parameters @f$\alpha@f$ and
    * @f$\beta@f$, over the interval @f$[a, b]@f$.
    *  @return the variance of the beta distribution @f$\mbox{Var}[X] =
    * \alpha\beta/ [(\alpha+ \beta)^2 (\alpha+ \beta+ 1)]@f$.
    */
   public static double getVariance (double alpha, double beta, double a,
                                     double b) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");

      return ((alpha * beta)*(b-a)*(b-a)) /
              ((alpha + beta) * (alpha + beta) * (alpha + beta + 1));
   }

   /**
    * Computes the standard deviation of the beta distribution with
    * parameters @f$\alpha@f$ and @f$\beta@f$, over the interval @f$[0,
    * 1]@f$.
    *  @return the standard deviation of the Beta distribution
    */
   public static double getStandardDeviation (double alpha, double beta) {
      return Math.sqrt (BetaDist.getVariance (alpha, beta));
   }

   /**
    * Computes the standard deviation of the beta distribution with
    * parameters @f$\alpha@f$ and @f$\beta@f$, over the interval @f$[a,
    * b]@f$.
    *  @return the standard deviation of the Beta distribution
    */
   public static double getStandardDeviation (double alpha, double beta,
                                              double a, double b) {
      return Math.sqrt (BetaDist.getVariance (alpha, beta, a, b));
   }

   /**
    * Returns the parameter @f$\alpha@f$ of this object.
    */
   public double getAlpha() {
      return alpha;
   }

   /**
    * Returns the parameter @f$\beta@f$ of this object.
    */
   public double getBeta() {
      return beta;
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

   @Deprecated
   public void setParams (double alpha, double beta,
                          double a, double b, int d) {
      setParams (alpha, beta, a, b);
   //   this.decPrec = d;
    }

/**
 * Sets the parameters of the current distribution. See the constructor.
 */
public void setParams (double alpha, double beta, double a, double b) {
      if (alpha <= 0.0)
         throw new IllegalArgumentException ("alpha <= 0");
      if (beta <= 0.0)
         throw new IllegalArgumentException ("beta <= 0");
      if (a >= b)
         throw new IllegalArgumentException ("a >= b");
      this.alpha = alpha;
      this.beta = beta;
      supportA = this.a = a;
      supportB = this.b = b;
      bminusa = b - a;
      double temp = Num.lnGamma (alpha);
      if (alpha == beta)
         temp *= 2.0;
      else
         temp += Num.lnGamma (beta);
      logBeta = temp - Num.lnGamma (alpha + beta);
      Beta = Math.exp(logBeta);
//      this.factor = 1.0 / (Beta * Math.pow (bminusa, alpha + beta - 1));
      this.logFactor = - logBeta - Math.log (bminusa) * (alpha + beta - 1);
    }

   /**
    * Return an array containing the parameters of the current
    * distribution as [@f$\alpha@f$, @f$\beta@f$, @f$a@f$, @f$b@f$].
    */
   @Override
   public double[] getParams () {
      double[] retour = {alpha, beta, a, b};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   @Override
   public String toString () {
      return getClass().getSimpleName() + " : alpha = " + alpha + ", beta = " + beta;
   }

}