/*
 * Class:        ChiSquareNoncentralGamGen
 * Description:  noncentral chi-square random variate generators using the
                 additive property of the noncentral chi-square distribution
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
 * using the additive property of the noncentral chi square distribution
 * @cite tKRI06a&thinsp;. It uses the following algorithm: generate a real
 * @f$X \sim N(\sqrt{\lambda}, 1)@f$ from a normal distribution with
 * variance 1, generate a real @f$Y \sim\Gamma((\nu- 1)/2, 1/2)@f$ from a
 * gamma distribution, then return @f$X^2 + Y@f$. Here @f$\nu@f$ is the
 * number of degrees of freedom and @f$\lambda@f$ is the noncentrality
 * parameter.
 *
 * To generate the normal variates, one uses the fast *acceptance-complement
 * ratio* method in @cite rHOR90a&thinsp; (see class
 * @ref umontreal.ssj.randvar.NormalACRGen ). To generate the gamma variates,
 * one uses acceptance-rejection for @f$\alpha<1@f$, and
 * acceptance-complement for @f$\alpha\ge1@f$, as proposed in
 * @cite rAHR72b, @cite rAHR82a&thinsp; (see class
 * @ref umontreal.ssj.randvar.GammaAcceptanceRejectionGen ).
 *
 * This noncentral chi square generator is faster than the generator
 * @ref umontreal.ssj.randvar.ChiSquareNoncentralPoisGen
 * on the next page of this guide
 * . For small @f$\lambda@f$, it is nearly twice as fast. As @f$\lambda@f$
 * increases, it is still faster but not as much.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class ChiSquareNoncentralGamGen extends ChiSquareNoncentralGen {
   private double racLam = -1.0;

   /**
    * Creates a noncentral chi square random variate generator with with
    * @f$\nu= @f$ `nu` degrees of freedom and noncentrality parameter
    * @f$\lambda= @f$ `lambda` using stream `stream`, as described above.
    */
   public ChiSquareNoncentralGamGen (RandomStream stream,
                                     double nu, double lambda) {
      super (stream, null);
      setParams (nu, lambda);
      racLam = Math.sqrt(lambda);
   }


   public double nextDouble() {
      return gamGen (stream, nu, racLam);
   }

/**
 * Generates a variate from the noncentral chi square distribution with
 * parameters @f$\nu= @f$&nbsp;`nu` and @f$\lambda= @f$&nbsp;`lambda` using
 * stream `stream`, as described above.
 */
public static double nextDouble (RandomStream stream,
                                    double nu, double lambda) {
      double racLam = Math.sqrt(lambda);
      return gamGen (stream, nu, racLam);
   }


//>>>>>>>>>>>>>>>>>>>>  P R I V A T E    M E T H O D S   <<<<<<<<<<<<<<<<<<<<

   private static double gamGen (RandomStream s, double nu, double racLam) {
      // racLam = sqrt(lambda)
      double x = NormalACRGen.nextDouble (s, racLam, 1.0);
      double y = GammaAcceptanceRejectionGen.nextDouble(s, 0.5*(nu - 1.0), 0.5);
      return x*x + y;
   }

}