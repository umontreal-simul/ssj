/*
 * Class:        JohnsonSLGen
 * Description:  random variate generators for the Johnson $S_L$ distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        July 2012

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