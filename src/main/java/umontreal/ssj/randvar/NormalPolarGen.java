/*
 * Class:        NormalPolarGen
 * Description:  normal random variate generators using the polar method
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
 * This class implements *normal* random variate generators using the *polar
 * method with rejection* @cite rMAR62a&thinsp;. Since the method generates
 * two variates at a time, the second variate is returned upon the next call
 * to  #nextDouble.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalPolarGen extends NormalGen {

   // used by polar method which calculate always two random values;
   private boolean available = false;
   private double[] variates = new double[2];
   private static double[] staticVariates = new double[2];

   /**
    * Creates a normal random variate generator with mean `mu` and
    * standard deviation `sigma`, using stream `s`.
    */
   public NormalPolarGen (RandomStream s, double mu, double sigma) {
      super (s, null);
      setParams (mu, sigma);
   }

   /**
    * Creates a standard normal random variate generator with @f$\mu=
    * 0@f$ and @f$\sigma=1@f$, using stream `s`.
    */
   public NormalPolarGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a random variate generator for the normal distribution
    * `dist` and stream `s`.
    */
   public NormalPolarGen (RandomStream s, NormalDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getSigma());
   }


   public double nextDouble() {
      if (available) {
         available = false;
         return mu + sigma*variates[1];
      }
      else {
         polar (stream, mu, sigma, variates);
         available = true;
         return mu + sigma*variates[0];
      }
   }

/**
 * Generates a variate from the normal distribution with parameters @f$\mu=
 * @f$&nbsp;`mu` and @f$\sigma= @f$&nbsp;`sigma`, using stream `s`.
 */
public static double nextDouble (RandomStream s, double mu, double sigma) {
      polar (s, mu, sigma, staticVariates);
      return mu + sigma*staticVariates[0];
   }
//>>>>>>>>>>>>>>>>>>>>  P R I V A T E     M E T H O D S   <<<<<<<<<<<<<<<<<<<<
   // Polar method with rejection
   private static void polar (RandomStream stream, double mu,
                              double sigma, double[] variates) {
      double x, y, s;
      do {
        x = 2*stream.nextDouble() - 1;
        y = 2*stream.nextDouble() - 1;
        s = x*x + y*y;
      } while (s > 1.0 || s == 0.0);

      double temp = Math.sqrt (-2.0*Math.log (s)/s);

      variates[0] = y*temp;
      variates[1] = x*temp;
   }

}