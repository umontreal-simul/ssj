/*
 * Class:        HypoExponentialGen
 * Description:  random variate generators for the hypoexponential distribution 
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        January 2011

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
 * This class implements random variate generators for the *hypoexponential*
 * distribution (see classes  @ref umontreal.ssj.probdist.HypoExponentialDist
 * and  @ref umontreal.ssj.probdist.HypoExponentialDistQuick in package
 * `probdist` for the definition).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class HypoExponentialGen extends RandomVariateGen {

   /**
    * Creates a hypoexponential random variate generator with rates
    * @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i =
    * 1,…,k@f$, using stream `stream`.
    */
   public HypoExponentialGen (RandomStream stream, double[] lambda) {
      super (stream, new HypoExponentialDist(lambda));
   }

   /**
    * Creates a new generator for the hypoexponential distribution `dist`
    * with stream `stream`.
    */
   public HypoExponentialGen (RandomStream stream, HypoExponentialDist dist) {
      super (stream, dist);
//      if (dist != null)
//         setParams (dist.getLambda());
   }

   /**
    * Uses inversion to generate a new hypoexponential variate with rates
    * @f$\lambda_i = @f$ <tt>lambda[</tt>@f$i-1@f$<tt>]</tt>, @f$i =
    * 1,…,k@f$, using stream `stream`. The inversion uses a root-finding
    * method and is very slow.
    */
   public static double nextDouble (RandomStream stream, double[] lambda) {
      return HypoExponentialDist.inverseF (lambda, stream.nextDouble());
   }

   /**
    * Returns the @f$\lambda_i@f$ associated with this object.
    */
   public double[] getLambda() {
      return ((HypoExponentialDist)dist).getLambda();
   }

   /**
    * Sets the rates @f$\lambda_i = @f$ <tt>lam[</tt>@f$i-1@f$<tt>]</tt>,
    * @f$i = 1,…,k@f$ of this object.
    */
   public void setLambda (double[] lambda) {
      ((HypoExponentialDist)dist).setLambda(lambda);
   }

}