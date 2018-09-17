/*
 * Class:        StudentPolarGen
 * Description:  Student-t random variate generators using the polar method
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
 * This class implements *Student* random variate generators using the
 * *polar* method of @cite rBAI94a&thinsp;. The code is adapted from UNURAN
 * (see @cite iLEY02a&thinsp;).
 *
 * The non-static `nextDouble` method generates two variates at a time and
 * the second one is saved for the next call. A pair of variates is generated
 * every second call. In the static case, two variates are generated per call
 * but only the first one is returned and the second is discarded.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class StudentPolarGen extends StudentGen {

    private boolean available = false;
    private double[] variates = new double[2];
    private static double[] staticVariates = new double[2];
    // Used by the polar method.

   /**
    * Creates a Student random variate generator with @f$n@f$ degrees of
    * freedom, using stream `s`.
    */
   public StudentPolarGen (RandomStream s, int n) {
      super (s, null);
      setN (n);
   }

   /**
    * Creates a new generator for the Student distribution `dist` and
    * stream `s`.
    */
   public StudentPolarGen (RandomStream s, StudentDist dist) {
      super (s, dist);
      if (dist != null)
         setN (dist.getN ());
   }
 

   public double nextDouble() {
       if (available) {
          available = false;
          return variates[1];
       }
       else {
          polar (stream, n, variates);
          available = true;
          return variates[0];
       }
   }

/**
 * Generates a new variate from the Student distribution with @f$n =
 * @f$&nbsp;`n` degrees of freedom, using stream `s`.
 */
public static double nextDouble (RandomStream s, int n) {
      polar (s, n, staticVariates);
      return staticVariates[0];
   }
//>>>>>>>>>>>>>>>>>>>>  P R I V A T E S    M E T H O D S   <<<<<<<<<<<<<<<<<<<<


/*****************************************************************************
 *                                                                           *
 * Student's t Distribution: Polar Method                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from Student's t distribution with  *
 *               parameters n > 0.                                          *
 *                                                                           *
 * REFERENCE:  - R.W. Bailey (1994): Polar generation of random variates     *
 *               with the t-distribution,                                    *
 *               Mathematics of Computation 62, 779-781.                     *
 *                                                                           *
 * Implemented by F. Niederl, 1994                                           *
 *****************************************************************************
 *                                                                           *
 * The polar method of Box/Muller for generating Normal variates is adapted  *
 * to the Student-t distribution. The two generated variates are not         *
 * independent and the expected no. of uniforms per variate is 2.5464.       *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 ****************************************************************************/

   private static void polar (RandomStream stream, int n, double[] variates) {
      double u,v,w;
      do {
         u = 2. * stream.nextDouble() - 1.;
         v = 2. * stream.nextDouble() - 1.;
         w = u*u + v*v;
       } while (w > 1.);

      double temp = Math.sqrt (n*(Math.exp (-2./n*Math.log (w)) - 1.)/w);
      variates[0] = u*temp;
      variates[1] = v*temp;
   }


}