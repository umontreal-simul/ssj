/*
 * Class:        BetaSymmetricalGen
 * Description:  random variate generators for the symmetrical beta distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
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