/*
 * Class:        FNoncentralGen
 * Description:  random variate generators for the noncentral F-distribution
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
 * This class implements random variate generators for the <em>noncentral
 * F</em>-distribution. If @f$X@f$ is a noncentral chi-square random variable
 * with @f$\nu_1 > 0@f$ degrees of freedom and noncentrality parameter
 * @f$\lambda> 0@f$, and @f$Y@f$ is a chi-square random variable
 * (statistically independent of @f$X@f$) with @f$\nu_2>0@f$ degrees of
 * freedom, then
 * @f[
 *   F’ = \frac{X/\nu_1}{Y/\nu_2}
 * @f]
 * has a noncentral @f$F@f$-distribution.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class FNoncentralGen extends RandomVariateGen {
   private ChiSquareNoncentralGen noncenchigen;
   private ChiSquareGen chigen;
   private double nu1;   // degrees of freedom of noncenchigen
   private int nu2;   // degrees of freedom of chigen

   public double nextDouble()  {
      double x = noncenchigen.nextDouble();
      double y = chigen.nextDouble();
      return (x * nu2) / (y * nu1);
   }

   /**
    * Creates a *noncentral-F* random variate generator using noncentral
    * chi-square generator `ncgen` and chi-square generator `cgen`.
    */
   public FNoncentralGen (ChiSquareNoncentralGen ncgen, ChiSquareGen cgen) {
      super (null, null);
      setChiSquareNoncentralGen (ncgen);
      setChiSquareGen (cgen);
   }

   /**
    * Sets the noncentral chi-square generator to `ncgen`.
    */
   public void setChiSquareNoncentralGen (ChiSquareNoncentralGen ncgen) {
      nu1 = ncgen.getNu();
      noncenchigen = ncgen;
   }

   /**
    * Sets the chi-square generator to `cgen`.
    */
   public void setChiSquareGen (ChiSquareGen cgen) {
      nu2 = cgen.getN();
      chigen = cgen;
   }

}