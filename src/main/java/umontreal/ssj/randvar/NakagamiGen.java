/*
 * Class:        NakagamiGen
 * Description:  random variate generators for the Nakagami distribution.
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
 * This class implements random variate generators for the *Nakagami*
 * distribution. See the definition in
 * @ref umontreal.ssj.probdist.NakagamiDist of package `probdist`.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NakagamiGen extends RandomVariateGen {
   // Distribution parameters
   protected double a;
   protected double lambda;
   protected double c;

   /**
    * Creates a new Nakagami generator with parameters @f$a=@f$ `a`,
    * @f$\lambda=@f$ `lambda` and @f$c =@f$ `c`, using stream `s`.
    */
   public NakagamiGen (RandomStream s, double a, double lambda, double c) {
      super (s, new NakagamiDist (a, lambda, c));
      setParams (a, lambda, c);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public NakagamiGen (RandomStream s, NakagamiDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getA(), dist.getLambda(), dist.getC());
   }

   /**
    * Generates a variate from the *Nakagami* distribution with parameters
    * @f$a=@f$ `a`, @f$\lambda=@f$ `lambda` and @f$c =@f$ `c`, using
    * stream `s`.
    *  @param s            the random stream
    *  @param a            the location parameter
    *  @param lambda       the scale parameter
    *  @param c            the shape parameter
    *  @return Generates a variate from the *Nakagami* distribution
    */
   public static double nextDouble (RandomStream s, double a, double lambda,
                                    double c) {
      return NakagamiDist.inverseF (a, lambda, c, s.nextDouble());
   }

   /**
    * Returns the location parameter @f$a@f$ of this object.
    *  @return the location parameter mu
    */
   public double getA() {
      return a;
   }

   /**
    * Returns the scale parameter @f$\lambda@f$ of this object.
    *  @return the scale parameter mu
    */
   public double getLambda() {
      return lambda;
   }

   /**
    * Returns the shape parameter @f$c@f$ of this object.
    *  @return the shape parameter mu
    */
   public double getC() {
      return c;
   }


   protected void setParams (double a, double lambda, double c) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (c <= 0.0)
         throw new IllegalArgumentException ("c <= 0");
      this.a = a;
      this.lambda = lambda;
      this.c = c;
   }
}