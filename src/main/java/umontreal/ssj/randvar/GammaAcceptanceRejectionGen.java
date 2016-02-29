/*
 * Class:        GammaAcceptanceRejectionGen
 * Description:  gamma random variate generators using a method that combines
                 acceptance-rejection with acceptance-complement
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
import cern.jet.stat.Gamma;

/**
 * This class implements *gamma* random variate generators using a method
 * that combines acceptance-rejection with acceptance-complement, and
 * proposed in @cite rAHR72b, @cite rAHR82a&thinsp;. It uses
 * acceptance-rejection for @f$\alpha<1@f$ and acceptance-complement for
 * @f$\alpha\ge1@f$.
 *
 * The above papers may lead to numerical errors when @f$\alpha@f$ is small.
 * When @f$\alpha< 0.1@f$, the acceptance-rejection algorithm of
 * @cite rLIU13a&thinsp; will be used. This algorithm generates the log value
 * of a gamma variate and it can handle a very small @f$\alpha@f$ parameter.
 * Use the `nextDoubleLog` methods to get the natural log value of the
 * generated gamma.
 *
 *  When @f$\alpha@f$ is close to 0, the generated gamma variate may be
 * numerically too small to be represented by the primitive date type
 * `double`. The method `nextDoubleLog` can be used to generate the natural
 * log value of a gamma variate. This allows the generation of very small
 * gamma variates.
 *
 * For each gamma variate, the first uniform required is taken from the main
 * stream and all additional uniforms (after the first rejection) are
 * obtained from the auxiliary stream.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class GammaAcceptanceRejectionGen extends GammaGen {
    
   private RandomStream auxStream;
   private static final double q1 =  0.0416666664;
   private static final double q2 =  0.0208333723;
   private static final double q3 =  0.0079849875;
   private static final double q4 =  0.0015746717;
   private static final double q5 = -0.0003349403;
   private static final double q6 =  0.0003340332;
   private static final double q7 =  0.0006053049;
   private static final double q8 = -0.0004701849;
   private static final double q9 =  0.0001710320;
   private static final double a1 =  0.333333333;
   private static final double a2 = -0.249999949;
   private static final double a3 =  0.199999867;
   private static final double a4 = -0.166677482;
   private static final double a5 =  0.142873973;
   private static final double a6 = -0.124385581;
   private static final double a7 =  0.110368310;
   private static final double a8 = -0.112750886;
   private static final double a9 =  0.104089866;
   private static final double e1 =  1.000000000;
   private static final double e2 =  0.499999994;
   private static final double e3 =  0.166666848;
   private static final double e4 =  0.041664508;
   private static final double e5 =  0.008345522;
   private static final double e6 =  0.001353826;
   private static final double e7 =  0.000247453;

   private int gen;
   // UNURAN parameters for the distribution
   private double beta;
   private double gamma;
   // Generator parameters
   // Acceptance rejection combined with acceptance complement
   private static final int gs = 0;
   private static final int gd = 1;
   private double b;
   private double ss;
   private double s;
   private double d;
   private double r;
   private double q0;
   private double c;
   private double si;

   // Threshold on shape parameter alpha.
   // Method of Liu, Martin and Syring, @cite rLIU13a will be used if alpha
   // is smaller than this threshold.
   private static final double USE_LMS_THRESHOLD = 0.1;

   /**
    * Creates a gamma random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\lambda=@f$ `lambda`, using main
    * stream `s` and auxiliary stream `aux`. The auxiliary stream is used
    * when a random number of uniforms is required for a rejection-type
    * generation method.
    */
   public GammaAcceptanceRejectionGen (RandomStream s, RandomStream aux,
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
   public GammaAcceptanceRejectionGen (RandomStream s,
                                       double alpha, double lambda) {
      this (s, s, alpha, lambda);
   }

   /**
    * Creates a new generator object for the gamma distribution `dist`,
    * using main stream `s` and auxiliary stream `aux`. The auxiliary
    * stream is used when a random number of uniforms is required for a
    * rejection-type generation method.
    */
   public GammaAcceptanceRejectionGen (RandomStream s, RandomStream aux, 
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
   public GammaAcceptanceRejectionGen (RandomStream s, GammaDist dist) {
      this (s, s, dist);
   }

   /**
    * Returns the auxiliary stream associated with this object.
    */
   public RandomStream getAuxStream() {
      return auxStream;
   }

   /**
    * Generates a new gamma variate with parameters @f$\alpha=
    * @f$&nbsp;`alpha` and @f$\lambda= @f$&nbsp;`lambda`, using main
    * stream `s` and auxiliary stream `aux`.
    */
   public static double nextDouble (RandomStream s, RandomStream aux, 
                                    double alpha, double lambda) {
      int gen = gs;
      double s_ = 0, ss = 0, d = 0, q0 = 0, r = 0, c = 0, si = 0, b = 0;

         if (alpha < USE_LMS_THRESHOLD) {
            // algorithm of Liu, Martin and Syring (2015)
            return acceptanceRejectionLms(s, aux, alpha, 1.0/lambda);
         }
         // Code taken from UNURAN
         else if (alpha < 1.0) {
            gen = gs;
            b = 1.0 + 0.36788794412*alpha;       // Step 1
         }
         else {
            gen = gd;
            // Step 1. Preparations
            ss = alpha - 0.5;
            s_ = Math.sqrt (ss);
            d = 5.656854249 - 12.0*s_;

            // Step 4. Set-up for hat case
            r = 1.0 / alpha;
            q0 = ((((((((q9*r + q8)*r + q7)*r + q6)*r + q5)*r + q4)*
                    r + q3)*r + q2)*r + q1)*r;
            if (alpha > 3.686) {
              if (alpha > 13.022) {
                b = 1.77;
                si = 0.75;
                c = 0.1515/s_;
              }
              else {
                b = 1.654 + 0.0076 * ss;
                si = 1.68/s_ + 0.275;
                c = 0.062/s_ + 0.024;
              }
            }
            else {
              b = 0.463 + s_ - 0.178*ss;
              si = 1.235;
              c = 0.195/s_ - 0.079 + 0.016*s_;
            }
         }
         return acceptanceRejection
             (s, aux, alpha, 1.0/lambda, 0, gen, b, s_, ss, d, r, q0, c, si);
   }
public double nextDouble() {
      if (alpha < USE_LMS_THRESHOLD)
         return acceptanceRejectionLms(stream, auxStream, alpha, beta);

      return acceptanceRejection
       (stream, auxStream, alpha, beta, gamma, gen, b, s, ss, d, r, q0, c, si);
   }

/**
 * Same as `nextDouble (s, s, alpha, lambda)`.
 */
public static double nextDouble (RandomStream s, double alpha, 
                                    double lambda) {
      return nextDouble (s, s, alpha, lambda);
   }

   /**
    * Returns the natural log value of a new gamma variate. When the shape
    * parameter @f$\alpha@f$ is close to 0, generating the log value of
    * the gamma variate is numerically more stable than generating
    * directly the gamma value with `nextDouble`.
    */
   public double nextDoubleLog() {
      if (alpha < USE_LMS_THRESHOLD)
         return acceptanceRejectionLmsLog(stream, auxStream, alpha, beta);

      return Math.log(acceptanceRejection
       (stream, auxStream, alpha, beta, gamma, gen, b, s, ss, d, r, q0, c, si));
   }

   /**
    * Returns the natural log value of a new gamma variate with parameters
    * @f$\alpha= @f$&nbsp;`alpha` and @f$\lambda= @f$&nbsp;`lambda`,
    * using main stream `s` and auxiliary stream `aux`. When the shape
    * parameter @f$\alpha@f$ is close to 0, generating the log value of
    * the gamma variate is numerically more stable than generating
    * directly the gamma value with `nextDouble`.
    */
   public static double nextDoubleLog (RandomStream s, RandomStream aux, 
                                    double alpha, double lambda) {
      int gen = gs;
      double s_ = 0, ss = 0, d = 0, q0 = 0, r = 0, c = 0, si = 0, b = 0;

      if (alpha < USE_LMS_THRESHOLD) {
         // algorithm of Liu, Martin and Syring (2015)
         return acceptanceRejectionLmsLog(s, aux, alpha, 1.0/lambda);
      }
      else {
         return Math.log(nextDouble(s, aux, alpha, lambda));
      }
   }

   /**
    * Same as `nextDoubleLog (s, s, alpha, lambda)`.
    */
   public static double nextDoubleLog (RandomStream s, double alpha, 
                                    double lambda) {
      return nextDoubleLog (s, s, alpha, lambda);
   }


   private static double acceptanceRejection
      (RandomStream stream, RandomStream auxStream,
       double alpha, double beta, double gamma, int gen,
       double b, double s, double ss,
       double d, double r, double q0, double c, double si) {
      // Code taken from UNURAN
      double X, p, U, E;
      double q,sign_U,t,v,w,x;
      switch (gen) {
      case gs:
         while (true) {
            p = b*stream.nextDouble();
            stream = auxStream;
            if (p <= 1.0) {                   // Step 2. Case gds <= 1
               X = Math.exp (Math.log (p)/alpha);
               if (Math.log (stream.nextDouble()) <= -X)
                  break;
            }
            else {                           // Step 3. Case gds > 1
               X = -Math.log ((b - p) / alpha);
               if ( Math.log (stream.nextDouble()) <= ((alpha - 1.0)*Math.log (X)))
                  break;
            }
         }
         break;
      case gd:
        do {

            // Step 2. Normal deviate
            t = NormalDist.inverseF01 (stream.nextDouble());
            stream = auxStream;
            x = s + 0.5*t;
            X = x*x;
            if (t >= 0.)
               break;         // Immediate acceptance

            // Step 3. Uniform random number
            U = stream.nextDouble();
            if (d*U <= t*t*t) 
               break;         // Squeeze acceptance

            // Step 5. Calculation of q
            if (x > 0.) {
               // Step 6.
               v = t/(s + s);
               if (Math.abs (v) > 0.25)
                  q = q0 - s*t + 0.25*t*t + (ss + ss)*Math.log (1. + v);
               else
                  q = q0 + 0.5*t*t*((((((((a9*v + a8)*v + a7)*v + a6)*
                                     v + a5)*v + a4)*v + a3)*v + a2)*v + a1)*v;
               // Step 7. Quotient acceptance
               if (Math.log (1. - U) <= q)
                  break;
            }

            // Step 8. Double exponential deviate t
            while (true) {
               do {
                  E = -Math.log (stream.nextDouble());
                  U = stream.nextDouble();
                  U = U + U - 1.;
                  sign_U = (U > 0) ? 1. : -1.;
                  t = b + (E*si)*sign_U;
               } while (t <= -0.71874483771719);   // Step 9. Rejection of t

               // Step 10. New q(t)
               v = t/(s + s);
               if (Math.abs (v) > 0.25)
                  q = q0 - s*t + 0.25*t*t + (ss + ss)*Math.log (1. + v);
               else
                  q = q0 + 0.5*t*t * ((((((((a9*v + a8)*v + a7)*v + a6)*
                                          v + a5)*v + a4)*v + a3)*v + a2)*v + a1)*v;

               // Step 11.
               if (q <= 0.)
                  continue; 

               if (q > 0.5)
                  w = Math.exp (q) - 1.;
               else
                  w = ((((((e7 * q + e6) * q + e5) * q + e4) * q + e3) * q + e2) *
                         q + e1) * q;

               // Step 12. Hat acceptance
               if ( c * U * sign_U <= w*Math.exp (E - 0.5*t*t)) {
                  x = s + 0.5 * t;
                  X = x * x;
                  break;
               }
            }
         } while (false);
         break;
      default: throw new IllegalStateException();
      }

      return gamma + beta*X;
   }


   /**
    * Generates a gamma variable. Use only when the shape is small (less than 0.2).
    * The method of @cite rLIU2015 for small shape parameter.
    */
   private static double acceptanceRejectionLms(RandomStream stream, RandomStream auxStream,
      double alpha, double scale) {

      return Math.exp(acceptanceRejectionLmsLog(stream, auxStream, alpha, scale));  // return the gamma variable
   }


   /**
    * Generates the natural log of a gamma variable. Use only when the shape is very small.
    * When the shape parameter is small, the value of the gamma variable can be 
    * numerically too small to be represented in a double.
    * 
    * The method of @cite rLIU2015 for small shape parameter.
    */
   private static double acceptanceRejectionLmsLog(RandomStream stream, RandomStream auxStream,
      double alpha, double scale) {

      double lambda = 1.0 / alpha - 1.0;
      double w = alpha / Num.EBASE / (1.0 - alpha);
      double r = 1.0 / (1.0 + w);
      double c = 1.0 / Gamma.gamma(alpha + 1.0);

      double u = stream.nextDouble();
      double z = 0;
      double Z = 0;
      while(true) {
         if (u <= r) {
            z = - Math.log(u / r);
         }
         else {
            z = Math.log(auxStream.nextDouble()) / lambda;
         }
         if (functionHLms(z, alpha, c) / functionEtaLms(z, alpha, lambda, w, c) > auxStream.nextDouble()) {
            Z = z;
            break;
         }
         u = auxStream.nextDouble();
      }
      return Math.log(scale) - (Z / alpha); // return scaled log gamma variable
   }

   // The function h in Liu, Martin and Syring (2015).
   private static double functionHLms(double z, double alpha, double c) {
      return c * Math.exp(-z - Math.exp(-z / alpha));
   }

   // The function eta in Liu, Martin and Syring (2015).
   private static double functionEtaLms(double z, double alpha, double lambda, double w, double c) {
      if (z >= 0)
         return c * Math.exp(-z);
      else
       return c * w * lambda * Math.exp(lambda * z);
   }


   private void init() {
      // Code taken from UNURAN
      if (alpha < 1.0) {
         gen = gs;
         b = 1.0 + 0.36788794412*alpha;       // Step 1
      }
      else {
         gen = gd;
         // Step 1. Preparations
         ss = alpha - 0.5;
         s = Math.sqrt (ss);
         d = 5.656854249 - 12.0*s;

         // Step 4. Set-up for hat case
         r = 1.0 / alpha;
         q0 = ((((((((q9*r + q8)*r + q7)*r + q6)*r + q5)*r + q4)*
                 r + q3)*r + q2)*r + q1)*r;
         if (alpha > 3.686) {
           if (alpha > 13.022) {
             b = 1.77;
             si = 0.75;
             c = 0.1515/s;
           }
           else {
             b = 1.654 + 0.0076 * ss;
             si = 1.68/s + 0.275;
             c = 0.062/s + 0.024;
           }
         }
         else {
           b = 0.463 + s - 0.178*ss;
           si = 1.235;
           c = 0.195/s - 0.079 + 0.016*s;
         }
      }
   }

}