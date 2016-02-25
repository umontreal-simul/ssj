/*
 * Class:        UniformGen
 * Description:  random variate generators for the uniform distribution over the reals
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
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements random variate generators for the (continuous)
 * *uniform* distribution over the interval @f$(a,b)@f$, where @f$a@f$ and
 * @f$b@f$ are real numbers with @f$a < b@f$. The density is
 * @anchor REF_randvar_UniformGen_eq_funiform
 * @f[
 *   f(x) = 1/(b - a) \qquad\mbox{ for }a\le x\le b. \tag{funiform}
 * @f]
 * The (non-static) `nextDouble` method simply calls `inverseF` on the
 * distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class UniformGen extends RandomVariateGen {
   private double a;
   private double b;

   /**
    * Creates a uniform random variate generator over the interval
    * @f$(@f$<tt>a</tt>, <tt>b</tt>@f$)@f$, using stream `s`.
    */
   public UniformGen (RandomStream s, double a, double b) {
      super (s, new UniformDist(a, b));
      setParams (a, b);
   }

   /**
    * Creates a uniform random variate generator over the interval @f$(0,
    * 1)@f$, using stream `s`.
    */
   public UniformGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a new generator for the uniform distribution `dist` and
    * stream `s`.
    */
   public UniformGen (RandomStream s, UniformDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getA(), dist.getB());
   }

   /**
    * Generates a uniform random variate over the interval
    * @f$(@f$<tt>a</tt>, <tt>b</tt>@f$)@f$ by inversion, using stream `s`.
    */
   static public double nextDouble (RandomStream s, double a, double b) {
      return UniformDist.inverseF (a, b, s.nextDouble());
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
    * Sets the value of the parameters @f$a@f$ and @f$b@f$ for this
    * object.
    */
   private void setParams (double a, double b) {
      if (b <= a)
         throw new IllegalArgumentException ("b <= a");
      this.a = a;
      this.b = b;
   }
}