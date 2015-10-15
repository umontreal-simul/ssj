/*
 * Class:        JohnsonSBGen
 * Description:  random variate generators for the Johnson $S_B$ distribution
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author
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
 * This class implements random variate generators for the <em>Johnson
 * @f$S_B@f$</em> distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class JohnsonSBGen extends JohnsonSystemG {

   /**
    * Creates a JohnsonSB random variate generator.
    */
   public JohnsonSBGen (RandomStream s, double gamma, double delta,
                        double xi, double lambda) {
      super (s, new JohnsonSBDist(gamma, delta, xi, lambda));
      setParams (gamma, delta, xi, lambda);
   }

   /**
    * Creates a new generator for the JohnsonSB distribution `dist`, using
    * stream `s`.
    */
   public JohnsonSBGen (RandomStream s, JohnsonSBDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getGamma(), dist.getDelta(), dist.getXi(),
                    dist.getLambda());
   }

   /**
    * Uses inversion to generate a new JohnsonSB variate, using stream
    * `s`.
    */
   public static double nextDouble (RandomStream s, double gamma,
                                    double delta, double xi, double lambda) {
      return JohnsonSBDist.inverseF (gamma, delta, xi, lambda,
                                        s.nextDouble());
   }
}