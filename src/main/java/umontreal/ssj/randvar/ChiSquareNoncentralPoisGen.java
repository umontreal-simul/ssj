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