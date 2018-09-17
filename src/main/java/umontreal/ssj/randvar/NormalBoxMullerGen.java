/*
 * Class:        NormalBoxMullerGen
 * Description:  normal random variate generators using the Box-Muller method
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
 * This class implements *normal* random variate generators using the
 * *Box-Muller* method from @cite rBOX58a&thinsp;. Since the method generates
 * two variates at a time, the second variate is returned upon the next call
 * to the  #nextDouble.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalBoxMullerGen extends NormalGen {
   private boolean available = false;
   private double[] variates = new double[2];
   private static double[] staticVariates = new double[2];
   // used by polar method which calculate always two random values;

   /**
    * Creates a normal random variate generator with mean `mu` and
    * standard deviation `sigma`, using stream `s`.
    */
   public NormalBoxMullerGen (RandomStream s, double mu, double sigma) {
      super (s, null);
      setParams (mu, sigma);
   }

   /**
    * Creates a standard normal random variate generator with mean `0` and
    * standard deviation `1`, using stream `s`.
    */
   public NormalBoxMullerGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a random variate generator for the normal distribution
    * `dist` and stream `s`.
    */
   public NormalBoxMullerGen (RandomStream s, NormalDist dist) {
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
         boxMuller (stream, mu, sigma, variates);
         available = true;
         return mu + sigma*variates[0];
      }
   }

/**
 * Generates a variate from the normal distribution with parameters @f$\mu=
 * @f$&nbsp;`mu` and @f$\sigma= @f$&nbsp;`sigma`, using stream `s`.
 */
public static double nextDouble (RandomStream s, double mu, double sigma) {
      boxMuller (s, mu, sigma, staticVariates);
      return mu + sigma*staticVariates[0];
   }
//>>>>>>>>>>>>>>>>>>>>  P R I V A T E S    M E T H O D S   <<<<<<<<<<<<<<<<<<<< 

   private static void boxMuller (RandomStream stream, double mu, 
                                  double sigma, double[] variates) {
      final double pi = Math.PI;
      double u,v,s;
    
      u = stream.nextDouble();
      v = stream.nextDouble();
      s = Math.sqrt (-2.0*Math.log (u));
      variates[1] = s*Math.sin (2*pi*v);
      variates[0] = s*Math.cos (2*pi*v);
   }

}