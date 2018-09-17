/*
 * Class:        GammaRejectionLoglogisticGen
 * Description:  gamma random variate generators using a rejection method
                 with log-logistic envelopes
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
 * This class implements *gamma* random variate generators using a rejection
 * method with loglogistic envelopes, from @cite rCHE77a&thinsp;. For each
 * gamma variate, the first two uniforms are taken from the main stream and
 * all additional uniforms (after the first rejection) are obtained from the
 * auxiliary stream.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class GammaRejectionLoglogisticGen extends GammaGen {
    
   private RandomStream auxStream;

   // UNURAN parameters for the distribution
   private double beta;
   private double gamma;
   // Generator parameters
   // Rejection with log-logistic envelopes
   private double aa;
   private double bb;
   private double cc;

   /**
    * Creates a gamma random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\lambda=@f$ `lambda`, using main
    * stream `s` and auxiliary stream `aux`. The auxiliary stream is used
    * when a random number of uniforms is required for a rejection-type
    * generation method.
    */
   public GammaRejectionLoglogisticGen (RandomStream s, RandomStream aux,
                                        double alpha, double lambda) {
      super (s, null);
      auxStream = aux;
      setParams (alpha, lambda);
      beta  = 1.0/lambda;
      gamma = 0.0;
      init ();
   }

   /**
    * Creates a gamma random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\lambda=@f$ `lambda`, using stream
    * `s`.
    */
   public GammaRejectionLoglogisticGen (RandomStream s,
                                        double alpha, double lambda) {
      this (s, s, alpha, lambda);
   }

   /**
    * Creates a new generator object for the gamma distribution `dist`,
    * using main stream `s` and auxiliary stream `aux`. The auxiliary
    * stream is used when a random number of uniforms is required for a
    * rejection-type generation method.
    */
   public GammaRejectionLoglogisticGen (RandomStream s, RandomStream aux, 
                                        GammaDist dist) {
      super (s, dist);
      auxStream = aux;
      if (dist != null)
         setParams (dist.getAlpha(), dist.getLambda());
      beta  = 1.0/dist.getLambda();
      gamma = 0.0;
      init ();
   }

   /**
    * Creates a new generator object for the gamma distribution `dist` and
    * stream `s` for both the main and auxiliary stream.
    */
   public GammaRejectionLoglogisticGen (RandomStream s, GammaDist dist) {
      this (s, s, dist);
   }

   /**
    * Returns the auxiliary stream associated with this object.
    */
   public RandomStream getAuxStream() {
      return auxStream;
   }


   public double nextDouble() {
      return rejectionLogLogistic 
                      (stream, auxStream, alpha, beta, gamma, aa, bb, cc);
   }

/**
 * Generates a new gamma variate with parameters @f$\alpha= @f$&nbsp;`alpha`
 * and @f$\lambda= @f$&nbsp;`lambda`, using main stream `s` and auxiliary
 * stream `aux`.
 */
public static double nextDouble (RandomStream s, RandomStream aux, 
                                    double alpha, double lambda) {
      double aa, bb, cc;

      // Code taken from UNURAN
      aa = (alpha > 1.0) ? Math.sqrt (alpha + alpha - 1.0) : alpha;
      bb = alpha - 1.386294361;
      cc = alpha + aa;
    
      return rejectionLogLogistic (s, aux, alpha, 1.0/lambda, 0.0, aa, bb, cc);
   }

   /**
    * Same as  {@link #nextDouble() nextDouble(s, s, alpha, lambda)}.
    */
   public static double nextDouble (RandomStream s, double alpha, 
                                    double lambda) {
      return nextDouble (s, s, alpha, lambda);
   }


   private static double rejectionLogLogistic
      (RandomStream stream, RandomStream auxStream,
       double alpha, double beta, double gamma,
        double aa, double bb, double cc) {
      // Code taken from UNURAN
      double X;
      double u1,u2,v,r,z;

      while (true) {
         u1 = stream.nextDouble();
         u2 = stream.nextDouble();
         stream = auxStream;
         v = Math.log (u1/(1.0 - u1))/aa;
         X = alpha*Math.exp (v);
         r = bb + cc*v - X;
         z = u1*u1*u2;
         if (r + 2.504077397 >= 4.5*z) break;
         if (r >= Math.log (z)) break;
      }

      return gamma + beta*X;
   }

   private void init() {
      // Code taken from UNURAN
      aa = (alpha > 1.0) ? Math.sqrt (alpha + alpha - 1.0) : alpha;
      bb = alpha - 1.386294361;
      cc = alpha + aa;
   }

}