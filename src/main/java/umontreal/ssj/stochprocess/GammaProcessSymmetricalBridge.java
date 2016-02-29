/*
 * Class:        GammaProcessSymmetricalBridge
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @authors      Pierre Tremblay and Jean-Sébastien Parent
 * @since        July 2003
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
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.util.*;

/**
 * This class differs from `GammaProcessBridge` only in that it requires the
 * number of interval of the path to be a power of 2 and of equal size. It is
 * then possible to generate the bridge process using a special
 * implementation of the beta random variate generator (using the
 * *symmetrical* beta distribution) that is much faster (HOW MUCH? QUANTIFY!)
 * than the general case. Note that when the method `setObservationTimes` is
 * called, the equality of the size of the time steps is verified. To allow
 * for differences due to floating point errors, time steps are considered to
 * be equal if their relative difference is less than @f$10^{-15}@f$. <div
 * class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class GammaProcessSymmetricalBridge extends GammaProcessBridge {
    protected BetaSymmetricalGen BSgen;

   /**
    * Constructs a new `GammaProcessSymmetricalBridge` with parameters
    * @f$\mu= \mathtt{mu}@f$, @f$\nu= \mathtt{nu}@f$ and initial value
    * @f$S(t_0) = \mathtt{s0}@f$. The random variables are created using
    * the  @ref umontreal.ssj.rng.RandomStream `stream`. Note that the
    * same  @ref umontreal.ssj.rng.RandomStream `stream` is used for the
    * @ref umontreal.ssj.randvar.GammaGen and for the
    * @ref umontreal.ssj.randvar.BetaSymmetricalGen inluded in this class.
    */
   public GammaProcessSymmetricalBridge (double s0, double mu, double nu,
                                         RandomStream stream) {
        this (s0, mu, nu, new GammaGen (stream, new GammaDist (1.0)),
              new BetaSymmetricalGen (stream, new BetaSymmetricalDist (1.0)));
    }

   /**
    * Constructs a new `GammaProcessSymmetricalBridge` with parameters
    * @f$\mu= \mathtt{mu}@f$, @f$\nu= \mathtt{nu}@f$ and initial value
    * @f$S(t_0) = \mathtt{s0}@f$. Note that the
    * @ref umontreal.ssj.rng.RandomStream included in the
    * @ref umontreal.ssj.randvar.BetaSymmetricalGen is sets to the one
    * included in the  @ref umontreal.ssj.randvar.GammaGen to avoid
    * confusion. This  @ref umontreal.ssj.rng.RandomStream is then used to
    * generate all the random variables.
    */
   public GammaProcessSymmetricalBridge (double s0, double mu, double nu,
                                         GammaGen Ggen,
                                         BetaSymmetricalGen BSgen) {
        super (s0, mu, nu, Ggen, BSgen);
        this.BSgen = BSgen;
        BSgen.setStream(Ggen.getStream());
    }


    public double nextObservation()  {
      double s;
      if (bridgeCounter == -1) {
         s = x0 + Ggen.nextDouble(stream, mu2dTOverNu, muOverNu);
         if (s <= x0)
              s = setLarger (x0);
         bridgeCounter = 0;
         observationIndex = d;
      } else {
         int j = bridgeCounter * 3;
         int oldIndexL = wIndexList[j];
         int newIndex = wIndexList[j + 1];
         int oldIndexR = wIndexList[j + 2];

         double y = BSgen.nextDouble(stream, bMu2dtOverNuL[newIndex]);

         s = path[oldIndexL] + (path[oldIndexR] - path[oldIndexL]) * y;
         if (s <= path[oldIndexL])
             s = setLarger (path, oldIndexL, oldIndexR);
         bridgeCounter++;
         observationIndex = newIndex;
      }
      observationCounter = bridgeCounter + 1;
      path[observationIndex] = s;
      return s;
    }

    public double nextObservation (double nextT) {
        double s;
        if (bridgeCounter == -1) {
            t[d] = nextT;
            mu2dTOverNu = mu2OverNu * (t[d] - t[0]);
            s = x0 + Ggen.nextDouble(stream, mu2dTOverNu, muOverNu);
            if (s <= x0)
               s = setLarger (x0);
            bridgeCounter    = 0;
            observationIndex = d;
        } else {
            int j = bridgeCounter*3;
            int oldIndexL = wIndexList[j];
            int newIndex  = wIndexList[j + 1];
            int oldIndexR = wIndexList[j + 2];

            t[newIndex] = nextT;
            bMu2dtOverNuL[newIndex] = mu2OverNu
                                      * (t[newIndex] - t[oldIndexL]);

            double y =  BSgen.nextDouble(stream, bMu2dtOverNuL[newIndex]);

            s = path[oldIndexL] + (path[oldIndexR] - path[oldIndexL]) * y;
            if (s <= path[oldIndexL])
                s = setLarger (path, oldIndexL, oldIndexR);
            bridgeCounter++;
            observationIndex = newIndex;
        }
        observationCounter = bridgeCounter + 1;
        path[observationIndex] = s;
        return s;
    }

    public double[] generatePath() {
        int oldIndexL, oldIndexR, newIndex;

        path[d] = x0 + Ggen.nextDouble(stream, mu2dTOverNu, muOverNu);
        for (int j = 0; j < 3*(d-1); j+=3) {
            oldIndexL   = wIndexList[j];
            newIndex    = wIndexList[j + 1];
            oldIndexR   = wIndexList[j + 2];

            double y =  BSgen.nextDouble(stream, bMu2dtOverNuL[newIndex]);

            path[newIndex] = path[oldIndexL] +
              (path[oldIndexR] - path[oldIndexL]) * y;
            if (path[newIndex] <= path[oldIndexL])
                setLarger (path, oldIndexL, newIndex, oldIndexR);
        }
        observationIndex   = d;
        observationCounter = d;
        return path;
    }

    public double[] generatePath (double[] uniform01) {
        int oldIndexL, oldIndexR, newIndex;

        path[d] = x0 + GammaDist.inverseF(mu2dTOverNu, muOverNu, 10, uniform01[0]);
        for (int j = 0; j < 3*(d-1); j+=3) {
            oldIndexL   = wIndexList[j];
            newIndex    = wIndexList[j + 1];
            oldIndexR   = wIndexList[j + 2];

            double y =  BetaSymmetricalDist.inverseF(bMu2dtOverNuL[newIndex], uniform01[1 + j/3]);

            path[newIndex] = path[oldIndexL] +
              (path[oldIndexR] - path[oldIndexL]) * y;
           if (path[newIndex] <= path[oldIndexL])
               setLarger (path, oldIndexL, newIndex, oldIndexR);
        }
        observationIndex   = d;
        observationCounter = d;
        return path;
    }

    protected void init () {
        super.init ();
        if (observationTimesSet) {

            /* Testing to make sure number of observations n = 2^k */
            int k = 0;
            int x = d;
            int y = 1;
            while (x>1) {
            x = x / 2;
            y = y * 2;
            k++;
            }
            if (y != d) throw new IllegalArgumentException
            ( "GammaProcessSymmetricalBridge:"
                +"Number 'n' of observation times is not a power of 2" );

            /* Testing that time intervals are equidistant */
            boolean equidistant = true;
            double macheps = 1.0e-13; // Num.DBL_EPSILON;
            double dt = t[1] - t[0];
            for (int i=1; i<d; i++) {
                if ((t[i+1] - t[i]) != dt) { // not equidistant
                    equidistant = false;
                    /* This compensates the fact that the dt's
                    may be different due to numerical idiosyncracies */
                    if (dt != 0.0)
                        if (Math.abs ((t[i+1] - t[i]) - dt) / dt <= macheps)
                            equidistant = true;
                }
            }
            if (!equidistant) throw new IllegalArgumentException
                        ( "GammaProcessSymmetricalBridge:"
                        +"Observation times of sample paths are not equidistant" );
        }
    }
}