/*
 * Class:        NormalACRGen
 * Description:  normal random variate generators using the
                 acceptance-complement ratio method
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
 * *acceptance-complement ratio* method @cite rHOR90a&thinsp;. For all the
 * methods, the code was taken from UNURAN @cite iLEY02a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalACRGen extends NormalGen {

   /**
    * Creates a normal random variate generator with mean `mu` and
    * standard deviation `sigma`, using stream `s`.
    */
   public NormalACRGen (RandomStream s, double mu, double sigma) {
      super (s, null);
      setParams (mu, sigma);
   }

   /**
    * Creates a standard normal random variate generator with mean `0` and
    * standard deviation `1`, using stream `s`.
    */
   public NormalACRGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a random variate generator for the normal distribution
    * `dist` and stream `s`.
    */
   public NormalACRGen (RandomStream s, NormalDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getSigma());
   }
 
  
   public double nextDouble() {
      return nextDouble (stream, mu, sigma);
   }

/**
 * Generates a variate from the normal distribution with parameters @f$\mu=
 * @f$&nbsp;`mu` and @f$\sigma= @f$&nbsp;`sigma`, using stream `s`.
 */
public static double nextDouble (RandomStream s, double mu, double sigma) {
/***************************************************************************** 
*       SAMPLING A RANDOM NUMBER FROM THE                            
*       STANDARD NORMAL DISTRIBUTION N(0,1).
*---------------------------------------------------------------------------
* GENERATION METHOD :  Acceptance-Complement Ratio 
*--------------------------------------------------------------------------- 
* REFERENCE : 
* - UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien
*         
* -  W. Hoermann and G. Derflinger (1990):                       
*     The ACR Method for generating normal random variables,       
*     OR Spektrum 12 (1990), 181-185.                             
*****************************************************************************/ 
      final double  c1   = 1.448242853;
      final double  c2   = 3.307147487;
      final double  c3   = 1.46754004;
      final double  d1   = 1.036467755;
      final double  d2   = 5.295844968;
      final double  d3   = 3.631288474;
      final double  hm   = 0.483941449;
      final double  zm   = 0.107981933;
      final double  hp   = 4.132731354;
      final double  zp   = 18.52161694;
      final double  phln = 0.4515827053;
      final double  hm1  = 0.516058551;
      final double  hp1  = 3.132731354;
      final double  hzm  = 0.375959516;
      final double  hzmp = 0.591923442;
      final double  as   = 0.8853395638;
      final double  bs   = 0.2452635696;
      final double  cs   = 0.2770276848;
      final double  b    = 0.5029324303;
      final double  x0   = 0.4571828819;
      final double  ym   = 0.187308492;
      final double  ss   = 0.7270572718; 
      final double  t    = 0.03895759111;
    
      double X;
      double rn,x,y,z;
    
      do {
         y = s.nextDouble();
         if (y > hm1) { X = hp*y - hp1; break; }
         else if (y < zm) {  
           rn = zp*y - 1;
           X = (rn > 0) ? (1 + rn) : (-1 + rn);
           break;
         } 
         else if (y < hm) {  
           rn = s.nextDouble();
           rn = rn - 1 + rn;
           z = (rn > 0) ? 2 - rn : -2 - rn;
           if ((c1 - y)*(c3 + Math.abs (z)) < c2) { X = z; break; }
           else {  
             x = rn*rn;
             if ((y + d1)*(d3 + x) < d2) { X = rn; break; }
             else if (hzmp - y < Math.exp (-(z*z + phln)/2)) { X = z; break; }
             else if (y + hzm < Math.exp (-(x + phln)/2)) { X = rn; break; }
           }
        }
      
        while (true) {
           x = s.nextDouble();
           y = ym*s.nextDouble();
           z = x0 - ss*x - y;
           if (z > 0)  rn = 2 + y/x;
           else {
             x = 1 - x;
             y = ym - y;
             rn = -(2 + y/x);
           }
           if ((y - as + x)*(cs + x) + bs < 0) { X = rn; break; }
           else if (y < x + t)
              if (rn*rn < 4*(b - Math.log (x))) { X = rn; break; }
        }
      
      } while (false);
        
      return mu + sigma*X;

   }
 }