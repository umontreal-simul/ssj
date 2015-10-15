/*
 * Class:        TriangularGen
 * Description:  random variate generators for the triangular distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements random variate generators for the *triangular*
 * distribution. Its density is
 * @anchor REF_randvar_TriangularGen_eq_ftrian
 * @f[
 *   f(x) = \left\{\begin{array}{ll}
 *    \frac{2(x-a)}{(b-a)(m-a)} 
 *    & 
 *    \mbox{for } a\le x\le m, 
 *    \\ 
 *    \frac{2(b-x)}{(b-a)(b-m)} 
 *    & 
 *    \mbox{ for } m\le x\le b, 
 *    \\ 
 *    0 
 *    & 
 *    \mbox{ elsewhere, } 
 *   \end{array}\right. \tag{ftrian}
 * @f]
 * where @f$a\le m\le b@f$ (see, e.g., @cite sLAW00a&thinsp;).
 *
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class TriangularGen extends RandomVariateGen {
   private double a;
   private double b;
   private double m;

   /**
    * Creates a triangular random variate generator over the interval
    * (<tt>a</tt>, <tt>b</tt>), with parameter `m`, using stream `s`.
    */
   public TriangularGen (RandomStream s, double a, double b, double m) {
      super (s, new TriangularDist(a, b, m));
      setParams (a, b, m);
   }

   /**
    * Creates a triangular random variate generator over the interval
    * @f$(0, 1)@f$, with parameter `m`, using stream `s`.
    */
   public TriangularGen (RandomStream s, double m) {
      this (s, 0.0, 1.0, m);
   }

   /**
    * Creates a new generator for the triangular distribution `dist` and
    * stream `s`.
    */
   public TriangularGen (RandomStream s, TriangularDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getA(), dist.getB(), dist.getM());
   }

   /**
    * Generates a new variate from the triangular distribution with
    * parameters @f$a = @f$&nbsp;`a`, @f$b = @f$&nbsp;`b` and @f$m =
    * @f$&nbsp;`m` and stream `s`, using inversion.
    */
   public static double nextDouble (RandomStream s, double a, 
                                    double b, double m) {
       // the code is taken and adapted from unuran
       // file /distributions/c_triangular_gen.c
       return TriangularDist.inverseF (a, b, m, s.nextDouble());
   }

   /**
    * Returns the value of @f$a@f$ for this object.
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the value of @f$b@f$ for this object.
    */
   public double getB() {
      return b;
   }

   /**
    * Returns the value of @f$m@f$ for this object.
    */
   public double getM() {
      return m;
   }

   /**
    * Sets the value of the parameters @f$a@f$, @f$b@f$ and @f$m@f$ for
    * this object.
    */
   private void setParams (double a, double b, double m) {
      if ((a == 0.0 && b == 1.0) && (m < 0 || m > 1))
         throw new IllegalArgumentException ("m is not in [0,1]");
      else if (a >= b)
         throw new IllegalArgumentException ("a >= b");
      else if (m < a || m > b) 
         throw new IllegalArgumentException ("m is not in [a,b]");
      this.a = a;
      this.b = b;
      this.m = m;
   }
}