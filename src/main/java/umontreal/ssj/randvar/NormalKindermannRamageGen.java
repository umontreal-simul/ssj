/*
 * Class:        NormalKindermannRamageGen
 * Description:  normal random variate generators using the Kindermann-Ramage method
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
 * *Kindermann-Ramage* method @cite rKIN76a&thinsp;. The code was taken from
 * UNURAN @cite iLEY02a&thinsp;. It includes the correction of the error in
 * the original *Kindermann-Ramage* method found by the authors in
 * @cite rTIR04a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalKindermannRamageGen extends NormalGen {

   /**
    * Creates a normal random variate generator with mean `mu` and
    * standard deviation `sigma`, using stream `s`.
    */
   public NormalKindermannRamageGen (RandomStream s,
                                     double mu, double sigma) {
      super (s, null);
      setParams (mu, sigma);
   }

   /**
    * Creates a standard normal random variate generator with mean `0` and
    * standard deviation `1`, using stream `s`.
    */
   public NormalKindermannRamageGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a random variate generator for the normal distribution
    * `dist` and stream `s`.
    */
   public NormalKindermannRamageGen (RandomStream s, NormalDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getSigma());
   }


   public double nextDouble() {
      return kindermanRamage (stream, mu, sigma);
   }

/**
 * Generates a variate from the normal distribution with parameters @f$\mu=
 * @f$&nbsp;`mu` and @f$\sigma= @f$&nbsp;`sigma`, using stream `s`.
 */
public static double nextDouble (RandomStream s, double mu, double sigma) {
      return kindermanRamage (s, mu, sigma);
   }
//>>>>>>>>>>>>>>>>>>>>  P R I V A T E S    M E T H O D S   <<<<<<<<<<<<<<<<<<<<

  private static double kindermanRamage (RandomStream stream, double mu, double sigma) {
    final double XI = 2.216035867166471;
    final double PIhochK = 0.3989422804;
    double U, V, W, X;
    double t, z;

    U = stream.nextDouble();
    if (U < 0.884070402298758) {
      V = stream.nextDouble();
      X = XI*(1.131131635444180*U + V - 1.);
    }

    else if (U >= 0.973310954173898) {
      do {
        V = stream.nextDouble();
        W = stream.nextDouble();
        if (W == 0.) { t=0.; continue; }
        t = XI*XI/2. - Math.log (W);
      } while ( (V*V*t) > (XI*XI/2.) );
      X = (U < 0.986655477086949) ? Math.pow (2*t,0.5) : -Math.pow (2*t,0.5);
    }

    else if (U >= 0.958720824790463) {
      do {
        V = stream.nextDouble();
        W = stream.nextDouble();
        z = V - W;
        t = XI - 0.630834801921960*Math.min (V, W);
      } while (Math.max (V, W) > 0.755591531667601 && 0.034240503750111*
               Math.abs (z) > (PIhochK*Math.exp (t*t/(-2.)) -
                              0.180025191068563*(XI - Math.abs (t))) );
      X = (z < 0) ? t : -t;
    }

    else if (U >= 0.911312780288703) {
      do {
        V = stream.nextDouble();
        W = stream.nextDouble();
        z = V - W;
        t = 0.479727404222441 + 1.105473661022070*Math.min (V, W);
      } while (Math.max (V, W) > 0.872834976671790 && 0.049264496373128*
               Math.abs (z) > (PIhochK*Math.exp (t*t/(-2))
                              - 0.180025191068563*(XI - Math.abs (t))) );
      X = (z < 0) ? t : -t;
    }

    else {
      do {
        V = stream.nextDouble();
        W = stream.nextDouble();
        z = V - W;
        t = 0.479727404222441 - 0.595507138015940*Math.min (V, W);
        if (t < 0.0)
           continue;
      } while (Math.max (V, W) > 0.805777924423817 && 0.053377549506886*
               Math.abs (z) > (PIhochK*Math.exp (t*t/(-2)) -
                              0.180025191068563*(XI - Math.abs (t))) );
      X = (z < 0) ? t : -t;
    }

    return mu + sigma*X;
   }
}