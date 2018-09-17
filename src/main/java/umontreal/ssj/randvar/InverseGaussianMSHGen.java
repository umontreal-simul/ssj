/*
 * Class:        InverseGaussianMSHGen
 * Description:  inverse gaussian random variate generators
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
 * This class implements *inverse gaussian* random variate generators using
 * the many-to-one transformation method of Michael, Schucany and Haas (MHS)
 * @cite rMIC76a, @cite rDEV06a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class InverseGaussianMSHGen extends InverseGaussianGen {
   private NormalGen genN;

   /**
    * Creates an *inverse gaussian* random variate generator with
    * parameters @f$\mu= @f$ `mu` and @f$\lambda= @f$ `lambda`, using
    * streams `s` and `sn`.
    */
   public InverseGaussianMSHGen (RandomStream s, NormalGen sn,
                                 double mu, double lambda) {
      super (s, null);
      setParams (mu, lambda, sn);
   }

   /**
    * Creates a new generator for the distribution `dist` using streams
    * `s` and `sn`.
    */
   public InverseGaussianMSHGen (RandomStream s, NormalGen sn,
                                 InverseGaussianDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getLambda(), sn);
   }

   /**
    * Generates a new variate from the *inverse gaussian* distribution
    * with parameters @f$\mu= @f$ `mu` and @f$\lambda= @f$ `lambda`,
    * using streams `s` and `sn`.
    */
   public static double nextDouble (RandomStream s, NormalGen sn,
                                    double mu, double lambda) {
      return mhs(s, sn, mu, lambda);
   }
 

   public double nextDouble() {
      return mhs(stream, genN, mu, lambda);
   }


//>>>>>>>>>>>>>>>>>>>>  P R I V A T E     M E T H O D S   <<<<<<<<<<<<<<<<<<<<

   private static double mhs (RandomStream stream, NormalGen genN,
                              double mu, double lambda) {
      if (lambda <= 0.0)
         throw new IllegalArgumentException ("lambda <= 0");
      if (mu <= 0.0)
         throw new IllegalArgumentException ("mu <= 0");

      double z = genN.nextDouble ();
      double ymu = mu*z*z;
      double x1 = mu + 0.5*mu*ymu/lambda - 0.5*mu/lambda *
                  Math.sqrt(4.0*ymu*lambda + ymu*ymu);
      double u = stream.nextDouble();
      if (u <= mu/(mu + x1))
         return x1;
      return mu*mu/x1;
   }


   protected void setParams (double mu, double lambda, NormalGen sn) {
      setParams (mu, lambda);
      this.genN = sn;
   }

}