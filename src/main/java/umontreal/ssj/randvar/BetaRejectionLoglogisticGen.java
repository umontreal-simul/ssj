/*
 * Class:        BetaRejectionLoglogisticGen
 * Description:  beta random variate generators using the rejection method
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
import umontreal.ssj.util.Num;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * Implements *Beta* random variate generators using the rejection method
 * with log-logistic envelopes from @cite rCHE78a&thinsp;. The method draws
 * the first two uniforms from the main stream and uses the auxiliary stream
 * for the remaining uniforms, when more than two are needed (i.e., when
 * rejection occurs).  The current implementation is adapted from UNURAN.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class BetaRejectionLoglogisticGen extends BetaGen {
    
   private RandomStream auxStream;
   // Parameters for rejection with log-logistic envelopes
   private static final int bb = 0;
   private static final int bc = 1;
   private double am;
   private double bm;
   private double al;
   private double alnam;
   private double be;
   private double ga;
   private double si;
   private double rk1;
   private double rk2;

   /**
    * Creates a beta random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta` over the interval
    * @f$(0,1)@f$, using main stream `s` and auxiliary stream `aux`. The
    * auxiliary stream is used when a random number of uniforms is
    * required for a rejection-type generation method.
    */
   public BetaRejectionLoglogisticGen (RandomStream s, RandomStream aux,
                                       double alpha, double beta) {
      super (s, null);
      auxStream = aux;
      setParams (alpha, beta, 0.0, 1.0);
      init();
   }

   /**
    * Creates a beta random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta`, over the interval
    * @f$(0,1)@f$, using stream `s`.
    */
   public BetaRejectionLoglogisticGen (RandomStream s,
                                       double alpha, double beta) {
      this (s, s, alpha, beta);
   }

   /**
    * Creates a beta random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta` over the interval
    * (<tt>a</tt>, <tt>b</tt>), using main stream `s` and auxiliary stream
    * `aux`. The auxiliary stream is used when a random number of uniforms
    * is required for a rejection-type generation method.
    */
   public BetaRejectionLoglogisticGen (RandomStream s, RandomStream aux,
          double alpha, double beta, double a, double b) {
      super (s, null);
      auxStream = aux;
      setParams (alpha, beta, a, b);
      init();
   }

   /**
    * Creates a beta random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta`, over the interval
    * (<tt>a</tt>, <tt>b</tt>), using stream `s`.
    */
   public BetaRejectionLoglogisticGen (RandomStream s,
          double alpha, double beta, double a, double b) {
      this (s, s, alpha, beta, a, b);
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s` and auxiliary stream `aux`. The main stream is used for the
    * first uniforms (before a rejection occurs) and the auxiliary stream
    * is used afterwards (after the first rejection).
    */
   public BetaRejectionLoglogisticGen (RandomStream s, RandomStream aux, 
                                       BetaDist dist) {
      super (s, dist);
      auxStream = aux;
      if (dist != null)
         setParams (dist.getAlpha(), dist.getBeta(), dist.getA(), dist.getB());
      init();
   }

   /**
    * Same as  {@link #BetaRejectionLoglogisticGen()
    * BetaRejectionLoglogisticGen(s, s, dist)}. The auxiliary stream used
    * will be the same as the main stream.
    */
   public BetaRejectionLoglogisticGen (RandomStream s, BetaDist dist) {
      this (s, s, dist);
   }

   /**
    * Returns the auxiliary stream associated with that object.
    */
   public RandomStream getAuxStream() {
      return auxStream;
   }


   private void init() {
      // Code taken from UNURAN
      if (p > 1.0 && q > 1.0) {
         gen = bb;
         am = (p < q) ? p : q;
         bm = (p > q) ? p : q;
         al = am + bm;
         be = Math.sqrt ((al - 2.0)/(2.0 * p * q - al));
         ga = am + 1.0 / be;
      }
      else {
         gen = bc;
         am = (p > q) ? p : q;
         bm = (p < q) ? p : q;
         al = am + bm;
         alnam = al * Math.log (al/am) - 1.386294361;
         be = 1.0 / bm;
         si = 1.0 + am - bm;
         rk1 = si * (0.013888889 + 0.041666667 * bm) / (am * be - 0.77777778);
         rk2 = 0.25 + (0.5 + 0.25 / si) * bm;
      }
   }

    
   public double nextDouble() {
     /************************************
      Previously: return rejectionLogLogistic();
      Now executes code directly (without calling anything)
      ***********************************/

      // The code was taken from UNURAN
      double X = 0.0;
      double u1,u2,v,w,y,z,r,s,t;
      RandomStream stream = this.stream;
      switch (gen) {
      case bb:
         /* -X- generator code -X- */
         while (true) {
            /* Step 1 */
            u1 = stream.nextDouble();
            u2 = stream.nextDouble();
            stream = auxStream;
            v = be*Math.log (u1/(1.0 - u1));
            w = am*Math.exp (v);
            z = u1*u1* u2;
            r = ga*v - 1.386294361;
            s = am + r - w;

            /* Step 2 */
            if (s + 2.609437912 < 5.0*z) {
               /* Step 3 */
               t = Math.log (z);
               if (s < t)
                  /* Step 4 */
                  if (r + al*Math.log (al/(bm + w)) < t) 
                     continue;
               }

            /* Step 5 */
            X = equalsDouble (am, p) ? w/(bm + w) : bm/(bm + w);
            break;
         }
         /* -X- end of generator code -X- */
         break;
      case bc:
         while (true) {
            /* Step 1 */
            u1 = stream.nextDouble();
            u2 = stream.nextDouble();
            stream = auxStream;

            if (u1 < 0.5) {
               /* Step 2 */
               y = u1*u2;
               z = u1*y;

               if ((0.25*u2 - y + z) >= rk1) 
                  continue;  /* goto 1 */

               /* Step 5 */
               v = be*Math.log (u1/(1.0 - u1));
               if (v > 80.0) {
                  if (alnam < Math.log (z))
                     continue;
                  X = equalsDouble (am, p) ? 1.0 : 0.0;
                  break;
               }
               else {
                  w = am*Math.exp (v);
                  if ((al*(Math.log (al/(bm + w)) + v) - 1.386294361) <
                                                                 Math.log (z))
                     continue;  /* goto 1 */

                  /* Step 6_a */
                  X = !equalsDouble (am, p) ? bm/(bm + w) : w/(bm + w);
                  break;
               }
            }
            else {
               /* Step 3 */
               z = u1*u1*u2;
               if (z < 0.25) {
                  /* Step 5 */
                  v = be*Math.log (u1/(1.0 - u1));
                  if (v > 80.0) {
                     X = equalsDouble (am, p) ? 1.0 : 0.0;
                     break;
                  }

                  w = am*Math.exp (v);
                  X = !equalsDouble (am, p) ? bm/(bm + w) : w/(bm + w);
                  break;
               }
               else {
                  if (z >= rk2)
                     continue;
                  v = be*Math.log (u1/(1.0 - u1));
                  if ( v > 80.0) {
                     if (alnam < Math.log (z))
                        continue;
                     X = equalsDouble (am, p) ? 1.0 : 0.0;
                     break;
                  }
                  w = am*Math.exp (v);
                  if ((al*(Math.log (al/(bm + w)) + v) - 1.386294361) < 
                                                                Math.log (z))
                     continue;  /* goto 1 */

                  /* Step 6_b */
                  X = !equalsDouble (am, p) ? bm/(bm + w) : w/(bm + w);
                  break;
               }
            }
         }
         break;
      default:
         throw new IllegalStateException();
      }

      return a + (b-a)*X;
   }

   private static boolean equalsDouble (double a, double b) {
      if (a == b)
         return true;
      double absa = Math.abs (a);
      double absb = Math.abs (b);
      return Math.abs (a - b) <= Math.min (absa, absb)*Num.DBL_EPSILON;
   }

}