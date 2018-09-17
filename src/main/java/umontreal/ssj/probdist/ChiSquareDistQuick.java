/*
 * Class:        ChiSquareDistQuick
 * Description:  chi-square distribution
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
 * Provides a variant of  @ref ChiSquareDist with faster but less accurate
 * methods. The non-static version of `inverseF` calls the static version.
 * This method is not very accurate for small @f$n@f$ but becomes better as
 * @f$n@f$ increases. The other methods are the same as in
 * @ref ChiSquareDist.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class ChiSquareDistQuick extends ChiSquareDist {

   /**
    * Constructs a chi-square distribution with `n` degrees of freedom.
    */
   public ChiSquareDistQuick (int n) {
      super (n);
   }


   public double inverseF (double u) {
      return inverseF (n, u);
   }

/**
 * Computes a quick-and-dirty approximation of @f$F^{-1}(u)@f$, where @f$F@f$
 * is the *chi-square* distribution with @f$n@f$ degrees of freedom. Uses the
 * approximation given in Figure L.24 of @cite sBRA87a&thinsp;  over most of
 * the range. For @f$u < 0.02@f$ or @f$u > 0.98@f$, it uses the approximation
 * given in @cite tGOL73a&thinsp;  for @f$n \ge10@f$, and returns `2.0 *`
 * {@link GammaDist.inverseF(double,int,double) inverseF(n/2, 6, u)} for @f$n
 * < 10@f$ in order to avoid the loss of precision of the above
 * approximations. When @f$n \ge10@f$ or @f$0.02 < u < 0.98@f$, it is
 * between 20 to 30 times faster than the same method in  @ref ChiSquareDist
 * for @f$n@f$ between @f$10@f$ and @f$1000@f$ and even faster for larger
 * @f$n@f$.
 *
 * Note that the number @f$d@f$ of decimal digits of precision generally
 * increases with @f$n@f$. For @f$n=3@f$, we only have @f$d = 3@f$ over most
 * of the range. For @f$n=10@f$, @f$d=5@f$ except far in the tails where @f$d
 * = 3@f$. For @f$n=100@f$, one has more than @f$d=7@f$ over most of the
 * range and for @f$n=1000@f$, at least @f$d=8@f$. The cases @f$n = 1@f$ and
 * @f$n = 2@f$ are exceptions, with precision of about @f$d=10@f$.
 */
public static double inverseF (int n, double u) {
      /*
       * Returns an approximation of the inverse of Chi square cdf
       * with n degrees of freedom.
       * As in Figure L.24 of P.Bratley, B.L.Fox, and L.E.Schrage.
       *         A Guide to Simulation Springer-Verlag,
       *         New York, second edition, 1987.
       */

      if (u < 0.0 || u > 1.0)
         throw new IllegalArgumentException ("u is not in [0,1]");
      if (u <= 0.0)
         return 0.0;
      if (u >= 1.0)
         return Double.POSITIVE_INFINITY;

      final double SQP5 = 0.70710678118654752440;
      final double DWARF = 0.1e-15;
      final double ULOW = 0.02;
      double Z, arg, v, ch, sqdf;

      if (n == 1) {
          Z = NormalDist.inverseF01 ((1.0 + u)/2.0);
          return Z*Z;

      } else if (n == 2) {
         arg = 1.0 - u;
         if (arg < DWARF)
            arg = DWARF;
         return -Math.log (arg)*2.0;

     } else if ((u > ULOW) && (u < 1.0 - ULOW)) {
        Z = NormalDist.inverseF01 (u);
        sqdf = Math.sqrt ((double)n);
        v = Z * Z;

        ch = -(((3753.0 * v + 4353.0) * v - 289517.0) * v -
           289717.0) * Z * SQP5 / 9185400;

        ch = ch / sqdf + (((12.0 * v - 243.0) * v - 923.0)
           * v + 1472.0) / 25515.0;

        ch = ch / sqdf + ((9.0 * v + 256.0) * v - 433.0)
           * Z * SQP5 / 4860;

        ch = ch / sqdf - ((6.0 * v + 14.0) * v - 32.0) / 405.0;
        ch = ch / sqdf + (v - 7.0) * Z * SQP5 / 9;
        ch = ch / sqdf + 2.0 * (v - 1.0) / 3.0;
        ch = ch / sqdf + Z / SQP5;
        return n * (ch / sqdf + 1.0);

     } else if (n >= 10) {
        Z = NormalDist.inverseF01 (u);
        v = Z * Z;
        double temp;
        temp = 1.0 / 3.0 + (-v + 3.0) / (162.0 * n) -
           (3.0 * v * v + 40.0 * v + 45.0) / (5832.0 * n * n) +
           (301.0 * v * v * v - 1519.0 * v * v - 32769.0 * v -
           79349.0) / (7873200.0 * n * n * n);
        temp *= Z * Math.sqrt (2.0 / n);

        ch = 1.0 - 2.0 / (9.0 * n) + (4.0 * v * v + 16.0 * v -
           28.0) / (1215.0 * n * n) + (8.0 * v * v * v + 720.0 * v * v +
           3216.0 * v + 2904.0) / (229635.0 * n * n * n) + temp;

        return n * ch * ch * ch;

    } else {
    // Note: this implementation is quite slow.
    // Since it is restricted to the tails, we could perhaps replace
    // this with some other approximation (series expansion ?)
        return 2.0*GammaDist.inverseF (n/2.0, 6, u);
    }
}

}