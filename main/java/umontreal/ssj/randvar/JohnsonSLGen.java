/*
 * Class:        JohnsonSLGen
 * Description:  random variate generators for the Johnson $S_L$ distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        July 2012
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
 * This class implements random variate generators for the <em>Johnson
 * @f$S_L@f$</em> distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class JohnsonSLGen extends JohnsonSystemG {

   /**
    * Creates a JohnsonSL random variate generator.
    */
   public JohnsonSLGen (RandomStream s, double gamma, double delta,
                        double xi, double lambda) {
      super (s, new JohnsonSLDist(gamma, delta, xi, lambda));
      setParams (gamma, delta, xi, lambda);
   }

   /**
    * Creates a new generator for the JohnsonSL distribution `dist`, using
    * stream `s`.
    */
   public JohnsonSLGen (RandomStream s, JohnsonSLDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getGamma(), dist.getDelta(), dist.getXi(),
                    dist.getLambda());
   }

   /**
    * Uses inversion to generate a new JohnsonSL variate, using stream
    * `s`.
    */
   public static double nextDouble (RandomStream s, double gamma,
                                    double delta, double xi, double lambda) {
      return JohnsonSLDist.inverseF (gamma, delta, xi, lambda,
                                        s.nextDouble());
   }
}