/*
 * Class:        BetaStratifiedRejectionGen
 * Description:  beta random variate generators using the stratified 
                 rejection/patchwork rejection method
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
 * This class implements *Beta* random variate generators using the
 * stratified rejection/patchwork rejection method from @cite rSAK83a,
 * @cite rSTA93a&thinsp;. This method draws one uniform from the main stream
 * and uses the auxiliary stream for any additional uniform variates that
 * might be needed.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class BetaStratifiedRejectionGen extends BetaGen {
    
   private RandomStream auxStream;
   private int gen;

   // Parameters for stratified rejection/patchwork rejection
   private static final int b00 = 2;
   private static final int b01 = 3;
   private static final int b01inv = 4;
   private static final int b1prs = 5;
   private double pint;
   private double qint;
   private double p_;
   private double q_;
   private double c;
   private double t;
   private double fp;
   private double fq;
   private double ml;
   private double mu;
   private double p1;
   private double p2;
   private double s;
   private double m;
   private double D;
   private double Dl;
   private double x1;
   private double x2;
   private double x4;
   private double x5;
   private double f1;
   private double f2;
   private double f4;
   private double f5;
   private double ll;
   private double lr;
   private double z2;
   private double z4;
   private double p3;
   private double p4;

   /**
    * Creates a beta random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta`, over the interval
    * @f$(0,1)@f$, using main stream `s` and auxiliary stream `aux`. The
    * auxiliary stream is used when a random number of uniforms is
    * required for a rejection-type generation method.
    */
   public BetaStratifiedRejectionGen (RandomStream s, RandomStream aux,
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
   public BetaStratifiedRejectionGen (RandomStream s,
                                       double alpha, double beta) {
      this (s, s, alpha, beta);
   }

   /**
    * Creates a beta random variate generator with parameters
    * @f$\alpha=@f$ `alpha` and @f$\beta=@f$ `beta`, over the interval
    * (<tt>a</tt>, <tt>b</tt>), using main stream `s` and auxiliary stream
    * `aux`. The auxiliary stream is used when a random number of uniforms
    * is required for a rejection-type generation method.
    */
   public BetaStratifiedRejectionGen (RandomStream s, RandomStream aux,
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
   public BetaStratifiedRejectionGen (RandomStream s,
          double alpha, double beta, double a, double b) {
      this (s, s, alpha, beta, a, b);
   }

   /**
    * Creates a new generator for the distribution `dist`, using the given
    * stream `s` and auxiliary stream `aux`. The auxiliary stream is used
    * when a random number of variates must be drawn from the main stream.
    */
   public BetaStratifiedRejectionGen (RandomStream s, RandomStream aux, 
                                      BetaDist dist) {
      super (s, dist);
      auxStream = aux;
      if (dist != null)
         setParams (dist.getAlpha(), dist.getBeta(), dist.getA(), dist.getB());
      init();
   }

   /**
    * Same as  {@link #BetaStratifiedRejectionGen()
    * BetaStratifiedRejectionGen(s, s, dist)}. The auxiliary stream used
    * will be the same as the main stream.
    */
   public BetaStratifiedRejectionGen (RandomStream s, BetaDist dist) {
      this (s, s, dist);
   }

   /**
    * Returns the auxiliary stream associated with this object.
    */
   public RandomStream getAuxStream() {
      return auxStream;
   }

    
   public double nextDouble() {
      /************************************
      Previously: return stratifiedRejection();
      Now executes code directly (without calling anything)
      ***********************************/

      // The code was taken from UNURAN
      double X = 0.0;
      double U, V, Z, W, Y;
      RandomStream stream = this.stream;
      switch (gen) {
      case b00:
         /* -X- generator code -X- */
         while (true) {
           U = stream.nextDouble()*p2;
           stream = auxStream;
            if (U <= p1) {                  /*  X < t       */
               Z = Math.exp (Math.log (U/p1)/p);
               X = t*Z;
               /* squeeze accept:   L(x) = 1 + (1 - q)x   */
               V = stream.nextDouble()*fq;
               if (V <= 1. - q_*X)
                  break;
               /* squeeze reject:   U(x) = 1 + ((1 - t)^(q-1) - 1)/t * x  */
               if (V <= 1. + (fq - 1.) * Z) {
                  /* quotient accept:  quot(x) = (1 - x)^(q-1) / fq       */
                  if (Math.log (V) <= q_*Math.log (1. - X))
                     break;
               }
            }
            else {          /*  X > t  */
               Z = Math.exp (Math.log ((U-p1)/(p2-p1) )/q);
               X = 1. - (1. - t)*Z;
               /* squeeze accept:   L(x) = 1 + (1 - p)(1 - x)            */
               V = stream.nextDouble()*fp;
               if (V <= 1.0 - p_*(1. - X))
                  break;
               /* squeeze reject: U(x) = 1 + (t^(p-1) - 1)/(1 - t) * (1 - x) */
               if (V <= 1.0 + (fp - 1.) * Z) {
                  /* quotient accept:  quot(x) = x^(p-1) / fp             */
                  if (Math.log (V) <= p_*Math.log (X))  
                     break;
               }
            }
         }
         /* -X- end of generator code -X- */
         break;
      case b01:
      case b01inv:
         /* -X- generator code -X- */
         while (true) {
           U = stream.nextDouble()*p2;
           stream = auxStream;
            if (U <= p1) {    /*  X < t                                 */
               Z = Math.exp (Math.log (U/p1)/pint);
               X = t*Z;
               /* squeeze accept:   L(x) = 1 + m1*x,  ml = -m1          */
               V = stream.nextDouble();
               if (V <= 1. - ml * X)
                  break;
               /* squeeze reject:   U(x) = 1 + m2*x,  mu = -m2 * t      */
               if (V <= 1. - mu * Z)
                  /* quotient accept:  quot(x) = (1 - x)^(q-1)          */
                 if (Math.log (V) <= q_*Math.log (1. - X))
                     break;
            }
            else {             /*  X > t                                */
               Z = Math.exp (Math.log ((U-p1)/(p2-p1)) / qint);
               X = 1. - (1. - t)*Z;
               /* squeeze accept:   L(x) = 1 + (1 - p)(1 - x)            */
               V = stream.nextDouble()*fp;
               if (V <= 1. - p_ * (1. - X))
                  break;
               /* squeeze reject: U(x) = 1 + (t^(p-1) - 1)/(1 - t) * (1 - x) */
               if (V <= 1. + (fp - 1.) * Z)
                  /* quotient accept:  quot(x) = (x)^(p-1) / fp          */
                  if (Math.log (V) <= p_*Math.log (X))
                     break;
            }
         }
         if (p>q)
            /* p and q has been swapped */
            X = 1. - X;
         /* -X- end of generator code -X- */
         break;
      case b1prs:
         while (true) {
           U = stream.nextDouble()*p4;
           stream = auxStream;
            if (U <= p1) {
               /* immediate accept:  x2 < X < m, - f(x2) < W < 0         */
               W = U/Dl - f2;
               if (W <= 0.0) {
                  X = m - U/f2;
                  break;
               }
               /* immediate accept:  x1 < X < x2, 0 < W < f(x1)          */
               if (W <= f1) {
                  X = x2 - W/f1*Dl;
                  break;
               }
               /* candidates for acceptance-rejection-test              */
               U = stream.nextDouble();
               V = Dl*U;
               X = x2 - V;
               Y = x2 + V;
               /* squeeze accept:    L(x) = f(x2) (x - z2) / (x2 - z2)   */
               if (W*(x2 - z2) <= f2*(X - z2))
                  break;
               V = f2 + f2 - W;
               if (V < 1.0) {
             /* squeeze accept: L(x) = f(x2) + (1 - f(x2))(x - x2)/(m - x2)  */
                  if (V <= f2 + (1. - f2)*U) {
                     X = Y;
                     break;
                  }
                  /* quotient accept:   x2 < Y < m,   W >= 2f2 - f(Y)     */
                  if (V <= Math.exp ( p_*Math.log (Y/m) 
                                    + q_*Math.log ((10. - Y)/(1.0 - m)) ) ) {
                     X = Y;
                     break;
                  }
               }
            }
            else if (U <= p2) {
               U -= p1;
               /* immediate accept:  m < X < x4, - f(x4) < W < 0  */
               W = U/D - f4;
               if (W <= 0.) {
                  X = m + U/f4;
                  break;
               }
               /* immediate accept:  x4 < X < x5, 0 < W < f(x5)    */
               if (W <= f5) {
                  X = x4 + W/f5 * D;
                  break;
               }
               /* candidates for acceptance-rejection-test     */
               U = stream.nextDouble();
               V = D*U;
               X = x4 + V;
               Y = x4 - V;
               /* squeeze accept:    L(x) = f(x4) (z4 - x) / (z4 - x4)  */
               if (W*(z4 - x4) <= f4*(z4 - X))
                  break;
               V = f4 + f4 - W;
               if (V < 1.0) {
              /* squeeze accept: L(x) = f(x4) + (1 - f(x4))(x4 - x)/(x4 - m) */
                  if (V <= f4 + (1.0 - f4)*U) {
                     X = Y;
                     break;
                  }
                  /* quotient accept:   m < Y < x4,   W >= 2f4 - f(Y)    */
                  if (V <= Math.exp ( p_*Math.log (Y/m) 
                                    + q_*Math.log ((1.0 - Y)/(1.0 - m)))) {
                     X = Y;
                     break;
                  }
               }
            }
            else if (U <= p3) {              /*      X < x1     */
               U = (U - p2)/(p3 - p2);
               Y = Math.log (U);
               X = x1 + ll*Y;
               if (X <= 0.0)                   /*      X > 0!!    */
                  continue; 
               W = U*stream.nextDouble();
               /* squeeze accept:  L(x) = f(x1) (x - z1) / (x1 - z1)   */
               /*      z1 = x1 - ll,   W <= 1 + (X - x1)/ll         */
               if (W <= 1.0 + Y)
                  break;
               W *= f1;
            }
            else {                      /*    x5 < X       */
               U = (U - p3)/(p4 - p3);
               Y = Math.log (U);
               X = x5 - lr*Y;
               if (X >= 1.0)                /*      X < 1!!    */
                  continue;
               W = U*stream.nextDouble();
               /* squeeze accept: L(x) = f(x5) (z5 - x) / (z5 - x5)    */
               /*                z5 = x5 + lr,   W <= 1 + (x5 - X)/lr  */
               if (W <= 1.0 + Y)
                  break;
               W *= f5;
            }
            /* density accept:  f(x) = (x/m)^(p_) ((1 - x)/(1 - m))^(q_) */
            if (Math.log (W) <= p_*Math.log (X/m) 
                                + q_*Math.log ((1.0 - X)/(1.0 - m)))
               break;
         }
         /* -X- end of generator code -X- */
         break;
      default:  throw new IllegalStateException();
      }

      return gen == b01inv ? a + (b-a)*(1.0 - X) : a + (b-a)*X;
   }

   public static double nextDouble (RandomStream s, 
                                    double alpha, double beta, 
                                    double a, double b) {
      return BetaDist.inverseF (alpha, beta, a, b, 15, s.nextDouble());
   }


   private void init() {
      // Code taken from UNURAN
      if (p > 1.) {
         if (q>1.)    /* p > 1 && q > 1 */
            gen = b1prs;
         else {        /* p > 1 && q <= 1 */
            gen = b01inv;
            double temp = p;
            p = q;
            q = temp;
         }
      }
      else {
         if (q>1.)    /* p <= 1 && q > 1 */
            gen = b01;
         else         /* p <= 1 && q <= 1 */
            gen = b00;
      }

      switch (gen) {
      case b00:
         /* -X- setup code -X- */
         p_ = p - 1.;
         q_ = q - 1.;
         c = (q*q_)/(p*p_);                              /* q(1-q) / p(1-p) */
         t = (c == 1.) ? 0.5 : (1. - Math.sqrt (c))/(1. - c);  /* t = t_opt */
         fp = Math.exp (p_*Math.log (t));
         fq = Math.exp (q_*Math.log (1. - t));      /* f(t) = fa * fb  */
  
         p1 = t/p;                                  /* 0 < X < t       */
         p2 = (1. - t)/q + p1;                    /* t < X < 1       */
         /* -X- end of setup code -X- */
         break;
      case b01:
      case b01inv:
         /* -X- setup code -X- */
         /* internal use of p and q */
         if (p > q) {
            /* swap p and q */
            pint = q;
            qint = p;
         }
         else {
            pint = p;
            qint = q;
         }

         p_ = pint - 1.;
         q_ = qint - 1.;
         t = p_/(pint - qint);        /* one step Newton * start value t   */
         fq = Math.exp ((q_ - 1.)*Math.log (1. - t));
         fp = pint - (pint + q_)*t;
         t -= (t - (1. - fp)*(1. - t)*fq/qint)/(1. - fp*fq);
         fp = Math.exp (p_*Math.log (t));
         fq = Math.exp (q_*Math.log (1. - t));     /* f(t) = fa * fb  */
         if (q_ <= 1.0) {
            ml = (1. - fq)/t;                      /*   ml = -m1      */
            mu = q_ * t;                          /*   mu = -m2 * t  */
         }
         else {
            ml = q_;
            mu = 1. - fq;
         }
         p1 = t/pint;                             /*  0 < X < t      */
         p2 = fq*(1. - t)/qint + p1;              /*  t < X < 1      */
         /* -X- end of setup code -X- */
         break;
      case b1prs:
         /* -X- setup code -X- */
         p_ = p - 1.0;
         q_ = q - 1.0;
         s = p_ + q_;
         m = p_/s;

         if (p_ > 1.0 || q_ > 1.0)
            D = Math.sqrt (m * (1. - m)/(s - 1.0));

         if (p_ <= 1.0) {
            x2 = (Dl = m * 0.5);
            x1 = z2 = f1 = ll = 0.0;
         }
         else {
            x2 = m - D;
            x1 = x2 - D;
            z2 = x2*(1.0 - (1.0 - x2)/(s*D));
            if (x1 <= 0.0 || (s - 6.0)*x2 - p_ + 3.0 > 0.0) {
               x1 = z2;  x2 = (x1 + m)*0.5;
               Dl = m - x2;
            }
            else {
               Dl = D;
            }
            f1 = Math.exp ( p_*Math.log (x1/m) 
                          + q_*Math.log ((1.0 - x1)/(1.0 - m)) );
            ll = x1*(1.0 - x1)/(s*(m - x1));            /* z1 = x1 - ll   */
         }
         f2 = Math.exp ( p_*Math.log (x2/m) 
                       + q_*Math.log ((1.0 - x2)/(1.0 - m)) );

         if (q_ <= 1.) {
            D = (1.0 - m)*0.5;
            x4 = 1.0 - D;
            x5 = z4 = 1.0;
            f5 = lr = 0.0;
         }
         else {
            x4 = m + D;
            x5 = x4 + D;
            z4 = x4*(1.0 + (1.0 - x4)/(s*D));
            if (x5 >= 1.0 || (s - 6.0)*x4 - p_ + 3.0 < 0.0) {
               x5 = z4;
               x4 = (m + x5)*0.5;
               D = x4 - m;
            }
            f5 = Math.exp ( p_*Math.log (x5/m) 
                          + q_*Math.log ((1.0 - x5)/(1. - m)) );
            lr = x5*(1.0 - x5)/(s*(x5 - m));            /* z5 = x5 + lr   */
         }
         f4 = Math.exp ( p_*Math.log (x4/m) 
                       + q_*Math.log ((1.0 - x4)/(1.0 - m)) );

         p1 = f2*(Dl + Dl);                                /*  x1 < X < m    */
         p2 = f4*(D  + D) + p1;                            /*  m  < X < x5   */
         p3 = f1*ll       + p2;                            /*       X < x1   */
         p4 = f5*lr       + p3;                            /*  x5 < X        */
         /* -X- end of setup code -X- */
         break;
      default: throw new IllegalStateException();
      }
   }

   private static boolean equalsDouble (double a, double b) {
      if (a == b)
         return true;
      double absa = Math.abs (a);
      double absb = Math.abs (b);
      return Math.abs (a - b) <= Math.min (absa, absb)*Num.DBL_EPSILON;
   }

}