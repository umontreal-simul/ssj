/*
 * Class:        HypergeometricGen
 * Description:  random variate generators for the hypergeometric distribution
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
 * This class implements random variate generators for the *hypergeometric*
 * distribution. Its mass function is
 *  (see, e.g., @cite rGEN98a&thinsp; (page 101))
 * @anchor REF_randvar_HypergeometricGen_eq_fheperg
 * @f[
 *   p(x) = \frac{ {m \choose x} {l - m\choose k-x}}{{l \choose k}} \qquad\mbox{for } x=\max(0,k-l+m), …, \min(k, m), \tag{fheperg}
 * @f]
 * where
 *  @f$m@f$, @f$l@f$ and @f$k@f$ are integers that satisfy @f$0< m\le l@f$
 * and @f$0 < k\le l@f$.
 *
 * The generation method is inversion using the chop-down algorithm
 * @cite sKAC85a&thinsp;
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class HypergeometricGen extends RandomVariateGenInt {
   private int m;
   private int l;
   private int k;

   /**
    * Creates a hypergeometric generator with parameters @f$m =
    * @f$&nbsp;`m`, @f$l = @f$&nbsp;`l` and @f$k = @f$&nbsp;`k`, using
    * stream `s`.
    */
   public HypergeometricGen (RandomStream s, int m, int l, int k) {
      super (s, new HypergeometricDist (m, l, k));
      setParams (m, l, k);
   }

   /**
    * Creates a new generator for distribution `dist`, using stream `s`.
    */
   public HypergeometricGen (RandomStream s, HypergeometricDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getM(), dist.getL(), dist.getK());
   }

   /**
    * Generates a new variate from the *hypergeometric* distribution with
    * parameters @f$m = @f$&nbsp;`m`, @f$l = @f$&nbsp;`l` and @f$k =
    * @f$&nbsp;`k`, using stream `s`.
    */
   public static int nextInt (RandomStream s, int m, int l, int k) {
      return HypergeometricDist.inverseF (m, l, k, s.nextDouble());
   }

   /**
    * Returns the @f$m@f$ associated with this object.
    */
   public int getM() {
      return m;
   }

   /**
    * Returns the @f$l@f$ associated with this object.
    */
   public int getL() {
      return l;
   }

   /**
    * Returns the @f$k@f$ associated with this object.
    */
   public int getK() {
      return k;
   }

   /**
    * Sets the parameter @f$n@f$ and @f$p@f$ of this object.
    */
   protected void setParams (int m, int l, int k) {
      if (l <= 0)
         throw new IllegalArgumentException ("l must be greater than 0");
      if (m <= 0 || m > l)
         throw new IllegalArgumentException ("m is invalid: 1<=m<l");
      if (k <= 0 || k > l)
         throw new IllegalArgumentException ("k is invalid: 1<=k<l");
      this.m = m;
      this.l = l;
      this.k = k;
   }
}