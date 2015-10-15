/*
 * Class:        ChiSquareNoncentralPoisGen
 * Description:  noncentral chi square random variate generators using Poisson
                 and central chi square generators
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
 * This class implements *noncentral chi square* random variate generators
 * using Poisson and central chi square generators. It uses the following
 * algorithm: generate a random integer @f$J
 * \sim\mbox{Poisson}(\lambda/2)@f$ from a Poisson distribution, generate a
 * random real @f$X \sim\Gamma(j + \nu/2, 1/2)@f$ from a gamma
 * distribution, then return @f$X@f$. Here @f$\nu@f$ is the number of
 * degrees of freedom and @f$\lambda@f$ is the noncentrality parameter.
 *
 * To generate the Poisson variates, one uses tabulated inversion for
 * @f$\lambda<10@f$, and the acceptance complement method for
 * @f$\lambda\ge10@f$, as in @cite rAHR82b&thinsp; (see class
 * @ref umontreal.ssj.randvar.PoissonTIACGen ). To generate the gamma
 * variates, one uses acceptance-rejection for @f$\alpha<1@f$, and
 * acceptance-complement for @f$\alpha\ge1@f$, as proposed in
 * @cite rAHR72b, @cite rAHR82a&thinsp; (see class
 * @ref umontreal.ssj.randvar.GammaAcceptanceRejectionGen ).
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ChiSquareNoncentralPoisGen extends ChiSquareNoncentralGen {

   /**
    * Creates a noncentral chi square random variate generator with
    * @f$\nu= @f$ `nu` degrees of freedom, and noncentrality parameter
    * @f$\lambda= @f$ `lambda`, using stream `stream` as described above.
    * Restriction: @f$\lambda\le4.29\times10^9@f$.
    */
   public ChiSquareNoncentralPoisGen (RandomStream stream,
                                      double nu, double lambda) {
      super (stream, null);
      setParams (nu, lambda);
      if (lambda > 4290000000.0)
         throw new UnsupportedOperationException("   lambda too large");
   }


   public double nextDouble() {
      return poisGenere (stream, nu, lambda);
   }

/**
 * Generates a variate from the noncentral chi square distribution with
 * parameters @f$\nu= @f$&nbsp;`nu`, and @f$\lambda= @f$&nbsp;`lambda`,
 * using stream `stream` as described above. Restriction:
 * @f$\lambda\le4.29\times10^9@f$.
 */
public static double nextDouble (RandomStream stream,
                                    double nu, double lambda) {
      return poisGenere (stream, nu, lambda);
   }


//>>>>>>>>>>>>>>>>>>>>  P R I V A T E S    M E T H O D S   <<<<<<<<<<<<<<<<<<<<

   private static double poisGenere (RandomStream s, double nu, double lambda) {
      int j = PoissonTIACGen.nextInt (s, 0.5*lambda);
      return GammaAcceptanceRejectionGen.nextDouble (s, 0.5*nu + j, 0.5);
   }

}