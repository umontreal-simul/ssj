/*
 * Class:        LogarithmicGen
 * Description:  random variate generators for the (discrete) logarithmic distribution
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
import umontreal.ssj.util.Num;

/**
 * This class implements random variate generators for the (discrete)
 * *logarithmic* distribution. Its mass function is
 * @anchor REF_randvar_LogarithmicGen_eq_flogar
 * @f[
 *   p(x) = \frac{-\theta^x}{x \log(1 - \theta)} \qquad\mbox{ for } x = 1,2,…, \tag{flogar}
 * @f]
 * where @f$0 < \theta<1@f$. It uses inversion with the LS chop-down
 * algorithm if @f$\theta< \theta_0@f$ and the LK transformation algorithm
 * if @f$\theta\ge\theta_0@f$, as described in @cite rKEM81a&thinsp;. The
 * threshold @f$\theta_0@f$ can be specified when invoking the constructor.
 * Its default value is @f$\theta_0 = 0.96@f$, as suggested in
 * @cite rKEM81a&thinsp;.
 *
 * A local copy of the parameter @f$\theta@f$ is maintained in this class.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_discrete
 */
public class LogarithmicGen extends RandomVariateGenInt {
   private static final double default_theta_limit = 0.96;

   private double theta_limit = default_theta_limit;
   private double theta;
   private double t;      // = log (1.0-theta).
   private double h;      // = -theta/log (1.0-theta).

   /**
    * Creates a logarithmic random variate generator with parameters
    * @f$\theta= @f$ `theta` and default value @f$\theta_0 = 0.96@f$,
    * using stream `s`.
    */
   public LogarithmicGen (RandomStream s, double theta) {
      this (s, theta, default_theta_limit);
   }

   /**
    * Creates a logarithmic random variate generator with parameters
    * @f$\theta= @f$ `theta` and @f$\theta_0 = \mathtt{theta0}@f$, using
    * stream `s`.
    */
   public LogarithmicGen (RandomStream s, double theta, double theta0) {
      super (s, null);
      this.theta = theta;
      theta_limit = theta0;
      init();
   }

   /**
    * Creates a new generator with distribution `dist` and stream `s`,
    * with default value @f$\theta_0 = 0.96@f$.
    */
   public LogarithmicGen (RandomStream s, LogarithmicDist dist) {
      this (s, dist, default_theta_limit);
   }

   /**
    * Creates a new generator with distribution `dist` and stream `s`,
    * with @f$\theta_0 = \mathtt{theta0}@f$.
    */
   public LogarithmicGen (RandomStream s, LogarithmicDist dist,
                          double theta0) {
      super (s, dist);
      theta_limit = theta0;
      if (dist != null)
         theta = dist.getTheta();
      init();
   }


   private void init () {
      if (theta <= 0.0 || theta >= 1.0)
         throw new IllegalArgumentException ("theta not in (0, 1)");
      if (theta >= theta_limit)
         h = Math.log1p(-theta);
      else
         t = -theta / Math.log1p(-theta);
   }

   public int nextInt() {
      if (theta < theta_limit)
         return ls (stream, theta, t);
      else   // Transformation
         return lk (stream, theta, h);
   }

   /**
    * Uses stream `s` to generate a new variate from the *logarithmic*
    * distribution with parameter @f$\theta=@f$ `theta`.
    */
   public static int nextInt (RandomStream s, double theta) {
      if (theta < default_theta_limit)
         return ls (s, theta, -theta/Math.log1p(-theta));
      else   // Transformation
         return lk (s, theta, Math.log1p(-theta));
   }



//>>>>>>>>>>>>>>>>>>>>  P R I V A T E    M E T H O D S   <<<<<<<<<<<<<<<<<<<<


   private static int ls (RandomStream stream, double theta, double t) {
      double u = stream.nextDouble();
      int x = 1;

      double p =  t;

      while (u > p) {
            u -= p;
            x++;
            p *= theta*((double) x - 1.0)/((double)x);
      }
      return x;
   }

   private static int lk (RandomStream stream, double theta, double h) {
      double u, v, p, q;
      int x;

      u = stream.nextDouble();
      if (u > theta)
            return 1;
      v = stream.nextDouble();
      q = 1.0 - Math.exp(v * h);
      if ( u <= q * q) {
           x = (int)(1. + (Math.log(u) / Math.log(q)));
           return x;
      }
      return ((u > q) ? 1 : 2);
   }

   /**
    * Returns the @f$\theta@f$ associated with this object.
    */
   public double getTheta() {
      return theta;
   }

   /**
    * Returns the @f$\theta_0@f$ associated with this object.
    */
   public double getTheta0() {
      return theta_limit;
   }

}