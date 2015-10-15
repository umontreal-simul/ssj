/*
 * Class:        StudentNoncentralGen
 * Description:  random variate generator for the noncentral Student-t distribution
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
 * This class implements random variate generators for the *noncentral
 * Student-t* distribution with @f$n>0@f$ degrees of freedom and
 * noncentrality parameter @f$\delta@f$. If @f$X@f$ is distributed according
 * to a normal distribution with mean @f$\delta@f$ and variance 1, and
 * @f$Y@f$ (statistically independent of @f$X@f$) is distributed according to
 * a chi-square distribution with @f$n@f$ degrees of freedom, then
 * @f[
 *   T’ = \frac{X}{\sqrt{Y/n}}
 * @f]
 * has a noncentral @f$t@f$-distribution with @f$n@f$ degrees of freedom and
 * noncentrality parameter @f$\delta@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class StudentNoncentralGen extends RandomVariateGen {
   private NormalGen normgen;
   private ChiSquareGen chigen;
   private int n;   // degrees of freedom of chi-square

   public double nextDouble()  {
      double x = normgen.nextDouble();
      double y = chigen.nextDouble();
      return x / Math.sqrt(y/n);
   }

   /**
    * Creates a *noncentral-t* random variate generator using normal
    * generator `ngen` and chi-square generator `cgen`.
    */
   public StudentNoncentralGen (NormalGen ngen, ChiSquareGen cgen) {
      super (null, null);
      setNormalGen (ngen);
      setChiSquareGen (cgen);
   }

   /**
    * Sets the normal generator to `ngen`.
    */
   public void setNormalGen (NormalGen ngen) {
      if (1.0 != ngen.getSigma())
         throw new IllegalArgumentException ("   variance of normal must be 1");
      normgen = ngen;
   }

   /**
    * Sets the chi-square generator to `cgen`.
    */
   public void setChiSquareGen (ChiSquareGen cgen) {
      chigen = cgen;
      n = cgen.getN();
   }

}