/*
 * Class:        ChiRatioOfUniformsGen
 * Description:  Chi random variate generators using the ratio of uniforms
                 method with shift
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
 * This class implements *Chi* random variate generators using the ratio of
 * uniforms method with shift.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ChiRatioOfUniformsGen extends ChiGen {

   /**
    * Creates a *chi* random variate generator with @f$\nu=@f$ `nu`
    * degrees of freedom, using stream `s`.
    */
   public ChiRatioOfUniformsGen (RandomStream s, int nu) {
      super (s, null);
      setParams (nu);
   }

   /**
    * Create a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public ChiRatioOfUniformsGen (RandomStream s, ChiDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getNu ());
   }

    
   public double nextDouble() {
      return ratioOfUniforms (stream, nu);
   }

   public static double nextDouble (RandomStream s, int nu) {
      if (nu <= 0)
         throw new IllegalArgumentException ("nu <= 0");
      return ratioOfUniforms (s, nu);
   }


//>>>>>>>>>>>>>>>>>>>>  P R I V A T E    M E T H O D S   <<<<<<<<<<<<<<<<<<<<
   
   //================================================================    
   // Method ratio of uniforms with shift.
   // code taken and adapted from unuran
   // file :  /distributions/c_chi_gen.c
   //=============================================================== 
   
   private static double ratioOfUniforms (RandomStream stream, int nu) {
      double u,v,z,zz,r;
      if (nu == 1) {
         while (true) {
           u = stream.nextDouble();
           v = stream.nextDouble() * 0.857763884960707;
           z = v / u;
           if (z < 0) continue;
           zz = z * z;
           r = 2.5 - zz;
           if (z < 0.)
              r = r + zz * z / (3. * z);
           if (u < r * 0.3894003915)
              break;
           if (zz > (1.036961043 / u + 1.4))
              continue;
           if (2 * Math.log(u) < (- zz * 0.5 ))
              break;
         }
      }

      else {  // nu > 1 
         final double b = Math.sqrt(nu - 1.);
         final double vm1 = - 0.6065306597 * (1. - 0.25 / (b * b + 1.));
         final double vm = (-b > vm1) ? -b : vm1;
         final double vp = 0.6065306597 * (0.7071067812 + b) / (0.5 + b);
         final double vd = vp - vm;
         while (true) {
           u = stream.nextDouble();
           v = stream.nextDouble() * vd + vm;
           z = v / u;
           if (z < -b)
              continue;
           zz = z * z;
           r = 2.5 - zz;
           if (z < 0.0)
           r = r + zz * z / (3.0 * (z + b));
           if (u < r * 0.3894003915) {
              z += b;
              break;
           }
           if (zz > (1.036961043 / u + 1.4))
              continue;
           if (2. * Math.log(u) < 
                      (Math.log(1.0 + z / b) * b * b - zz * 0.5 - z * b)) {
              z += b;
              break;
           }
         }
      } 
      return z;
   }
}