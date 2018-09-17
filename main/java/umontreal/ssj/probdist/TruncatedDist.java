/*
 * Class:        TruncatedDist
 * Description:  an arbitrary continuous distribution truncated
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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
import umontreal.ssj.functions.MathFunctionUtil;
import umontreal.ssj.functions.MathFunction;

/**
 * This container class takes an arbitrary continuous distribution and
 * truncates it to an interval @f$[a,b]@f$, where @f$a@f$ and @f$b@f$ can be
 * finite or infinite. If the original density and distribution function are
 * @f$f_0@f$ and @f$F_0@f$, the new ones are @f$f@f$ and @f$F@f$, defined by
 * @f[
 *   f(x) = \frac{f_0(x)}{F_0(b) - F_0(a)} \qquad\mbox{ for } a\le x\le b
 * @f]
 * and @f$f(x)=0@f$ elsewhere, and
 * @f[
 *   F(x) = \frac{F_0(x)-F_0(a)}{F_0(b)-F_0(a)} \qquad\mbox{ for } a\le x\le b.
 * @f]
 * The inverse distribution function of the truncated distribution is
 * @f[
 *   {F^{-1}}(u) = F_0^{-1}(F_0(a) + (F_0(b) - F_0(a))u)
 * @f]
 * where @f$F_0^{-1}@f$ is the inverse distribution function of the original
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_continuous
 */
public class TruncatedDist extends ContinuousDistribution {
   public static int NUMINTERVALS = 500;

   private ContinuousDistribution dist0;  // The original (non-truncated) dist.
   private double fa;                    // F(a)
   private double fb;                    // F(b)
   private double barfb;                 // bar F(b)
   private double fbfa;                  // F(b) - F(a)
   private double a;
   private double b;
   private double approxMean;
   private double approxVariance;
   private double approxStandardDeviation;

   /**
    * Constructs a new distribution by truncating distribution `dist` to
    * the interval @f$[a,b]@f$. Restrictions: @f$a@f$ and @f$b@f$ must be
    * finite.
    */
   public TruncatedDist (ContinuousDistribution dist, double a, double b) {
      setParams (dist, a, b);
   }
   public double density (double x) {
      if ((x < a) || (x > b))
         return 0;
      return dist0.density (x) / fbfa;
   }

   public double cdf (double x) {
      if (x <= a)
         return 0;
      else if (x >= b)
         return 1;
      else
         return (dist0.cdf (x) - fa) / fbfa;
   }

   public double barF (double x) {
      if (x <= a)
         return 1;
      else if (x >= b)
         return 0;
      else
         return (dist0.barF (x) - barfb) / fbfa;
   }

   public double inverseF (double u) {
      if (u == 0)
         return a;
      if (u == 1)
         return b;
      return dist0.inverseF (fa + fbfa * u);
   }

   /**
    * Returns an approximation of the mean computed with the Simpson
    * @f$1/3@f$ numerical integration rule.
    *  @exception UnsupportedOperationException the mean of the truncated
    * distribution is unknown
    */
   public double getMean() {
      if (Double.isNaN (approxMean))
         throw new UnsupportedOperationException("Undefined mean");
      return approxMean;
   }

   /**
    * Returns an approximation of the variance computed with the Simpson
    * @f$1/3@f$ numerical integration rule.
    *  @exception UnsupportedOperationException the mean of the truncated
    * distribution is unknown
    */
   public double getVariance() {
      if (Double.isNaN (approxVariance))
         throw new UnsupportedOperationException("Unknown variance");
      return approxVariance;
   }

   /**
    * Returns the square root of the approximate variance.
    *  @exception UnsupportedOperationException the mean of the truncated
    * distribution is unknown
    */
   public double getStandardDeviation() {
      if (Double.isNaN (approxStandardDeviation))
         throw new UnsupportedOperationException("Unknown standard deviation");
      return approxStandardDeviation;
   }

   /**
    * Returns the value of @f$a@f$.
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the value of @f$b@f$.
    */
   public double getB() {
      return b;
   }

   /**
    * Returns the value of @f$F_0(a)@f$.
    */
   public double getFa() {
      return fa;
   }

   /**
    * Returns the value of @f$F_0(b)@f$.
    */
   public double getFb() {
      return fb;
   }

   /**
    * Returns the value of @f$F_0(b) - F_0(a)@f$, the area under the
    * truncated density function.
    */
   public double getArea() {
      return fbfa;
   }

   /**
    * Sets the parameters `dist`, @f$a@f$ and @f$b@f$ for this object. See
    * the constructor for details.
    */
   public void setParams (ContinuousDistribution dist, double a, double b) {
      if (a >= b)
         throw new IllegalArgumentException ("a must be smaller than b.");
      this.dist0 = dist;
      if (a < dist.getXinf())
         a = dist.getXinf();
      if (b > dist.getXsup())
         b = dist.getXsup();
      supportA = this.a = a;
      supportB = this.b = b;
      fa = dist.cdf (a);
      fb = dist.cdf (b);
      fbfa = fb - fa;
      barfb = dist.barF (b);

      if (((a <= dist.getXinf()) && (b >= dist.getXsup())) ||
       ((a == Double.NEGATIVE_INFINITY) && (b == Double.POSITIVE_INFINITY))) {
         approxMean = dist.getMean();
         approxVariance = dist.getVariance();
         approxStandardDeviation = dist.getStandardDeviation();

      } else {
         // The mean is the integral of xf*(x) over [a, b].
         MomentFunction func1 = new MomentFunction (dist, 1);
         approxMean = MathFunctionUtil.simpsonIntegral (func1, a, b, NUMINTERVALS) / fbfa;

         // Estimate the integral of (x-E[X])^2 f*(x) over [a, b]
         MomentFunction func2 = new MomentFunction (dist, 2, approxMean);
         approxVariance = MathFunctionUtil.simpsonIntegral (func2, a, b, NUMINTERVALS) / fbfa;

         approxStandardDeviation = Math.sqrt (approxVariance);
      }
   }

   /**
    * Return a table containing the parameters of the current
    * distribution. This table is put in order: [@f$a@f$, @f$b@f$,
    * @f$F_0(a)@f$, @f$F_0(b)@f$, @f$F_0(b) - F_0(a)@f$].
    */
   public double[] getParams () {
      double[] retour = {a, b, fa, fb, fbfa};
      return retour;
   }

   /**
    * Returns a `String` containing information about the current
    * distribution.
    */
   public String toString () {
      return getClass().getSimpleName() + " : a = " + a + ", b = " + b + ", F(a) = " + fa + ", F(b) = " + fb + ", F(b)-F(a) = " + fbfa;
   }


   private static class MomentFunction implements MathFunction {
      private ContinuousDistribution dist;
      private int moment;
      private double offset;

      public MomentFunction (ContinuousDistribution dist, int moment) {
         this.dist = dist;
         this.moment = moment;
         offset = 0;
      }

      public MomentFunction (ContinuousDistribution dist, int moment, double offset) {
         this (dist, moment);
         this.offset = offset;
      }

      public double evaluate (double x) {
         double res = dist.density (x);
         final double offsetX = x - offset;
         for (int i = 0; i < moment; i++)
            res *= offsetX;
         return res;
      }
   }
}