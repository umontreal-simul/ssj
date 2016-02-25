/*
 * Class:        FisherFGen
 * Description:  random variate generators for the Fisher F distribution
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
 * This class implements random variate generators for the *Fisher F*
 * distribution with @f$n@f$ and @f$m@f$ degrees of freedom, where @f$n@f$
 * and @f$m@f$ are positive integers. The density function of this
 * distribution is
 * @anchor REF_randvar_FisherFGen_eq_FisherF
 * @f[
 *   f(x) = \frac{\Gamma({(n + m)}/2)n^{n/2}m^{m/2}}{\Gamma(n/2)\Gamma(m/2)} \frac{x^{{(n - 2)}/2}}{(m + nx)^{{(n + m)}/2}}, \qquad\mbox{for } x > 0 \tag{FisherF}
 * @f]
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class FisherFGen extends RandomVariateGen {
   protected int    n = -1;
   protected int    m = -1;

   /**
    * Creates a *Fisher F* random variate generator with @f$n@f$ and
    * @f$m@f$ degrees of freedom, using stream `s`.
    */
   public FisherFGen (RandomStream s, int n, int m) {
      super (s, new FisherFDist(n, m));
      setParams (n, m);
      }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public FisherFGen (RandomStream s, FisherFDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getN(), dist.getM());
   }

   /**
    * Generates a variate from the *Fisher F* distribution with @f$n@f$
    * and @f$m@f$ degrees of freedom, using stream `s`.
    */
   public static double nextDouble (RandomStream s, int n, int m) {
      return FisherFDist.inverseF (n, m, 15, s.nextDouble());
   }

   /**
    * Returns the parameter @f$n@f$ of this object.
    */
   public int getN() {
      return n;
   }

   /**
    * Returns the parameter @f$p@f$ of this object.
    */
   public int getM() {
      return m;
   }

   /**
    * Sets the parameters @f$n@f$ and @f$m@f$ of this object.
    */
   protected void setParams (int n, int m) {
      if (m <= 0)
         throw new IllegalArgumentException ("m <= 0");
      if (n <= 0)
         throw new IllegalArgumentException ("n <= 0");
      this.m = m;
      this.n = n;
   }
}