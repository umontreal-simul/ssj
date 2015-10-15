/*
 * Class:        BetaSymmetricalGen
 * Description:  random variate generators for the symmetrical beta distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */
package umontreal.ssj.randvar;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements random variate generators with the *symmetrical
 * beta* distribution with shape parameters @f$\alpha= \beta@f$, over the
 * interval @f$(0,1)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class BetaSymmetricalGen extends BetaGen {

   /**
    * Creates a new symmetrical beta generator with parameters
    * @f$\alpha=@f$ `alpha`, over the interval @f$(0,1)@f$, using stream
    * `s`.
    */
   public BetaSymmetricalGen (RandomStream s, double alpha) {
      this (s, new BetaSymmetricalDist (alpha));
   }

   /**
    * Creates a new generator for the distribution `dist`, using stream
    * `s`.
    */
   public BetaSymmetricalGen (RandomStream s, BetaSymmetricalDist dist) {
      super (s, dist);
   }
   public static double nextDouble (RandomStream s, double alpha) {
      return BetaSymmetricalDist.inverseF (alpha, s.nextDouble());
   }

}