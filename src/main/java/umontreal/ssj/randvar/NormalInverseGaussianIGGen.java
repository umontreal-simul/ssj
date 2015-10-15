/*
 * Class:        NormalInverseGaussianIGGen
 * Description:  normal inverse gaussian random variate generator
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        June 2008

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
 * This class implements a *normal inverse gaussian* (@f${NIG}@f$) random
 * variate generator by using a normal generator (@f$N@f$) and an inverse
 * gaussian generator (@f$IG@f$), as described in the following
 * @cite fWEB03a, @cite fKAL07a&thinsp;
 * @anchor REF_randvar_NormalInverseGaussianIGGen_nig2
 * @f{align}{
 *    Y 
 *    & 
 *   \sim
 *    {IG}(\delta/\gamma, \delta^2) \tag{nig2} 
 *    \\ 
 *   X \mid(Y=y) 
 *    & 
 *   \sim
 *    N(\mu+ \beta y, y). \nonumber
 * @f}
 * The normal @f$N(\mu, \sigma^2)@f$ has mean @f$\mu@f$ and variance
 * @f$\sigma^2@f$, while the inverse gaussian has the parametrization
 * described in
 * equation (
 * {@link REF_randvar_InverseGaussianGen_eq_fInverseGaussian
 * fInverseGaussian} ).
 *  If @f$\gamma= \sqrt{\alpha^2 - \beta^2}@f$ with @f$0 \le|\beta| <
 * \alpha@f$ and @f$\delta>0@f$, then @f$X \sim{NIG}(\alpha, \beta,
 * \mu, \delta)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalInverseGaussianIGGen extends NormalInverseGaussianGen {
   private NormalGen genN;
   private InverseGaussianGen genIG;

   /**
    * Creates a *normal inverse gaussian* random variate generator with
    * parameters @f$\alpha@f$, @f$\beta@f$ = `beta`, @f$\mu@f$ = `mu`
    * and @f$\delta@f$, using generators `ig` and `ng`, as described
    *  in eq. (
    * {@link REF_randvar_NormalInverseGaussianIGGen_nig2
    * nig2} ).
    *  The parameters @f$\alpha@f$ and @f$\delta@f$ are included in
    * generator `ig`.
    */
   public NormalInverseGaussianIGGen (InverseGaussianGen ig, NormalGen ng,
                                      double beta, double mu) {
      super (null, null);
      setParams (ig, ng, beta, mu);
   }

   /**
    * Generates a new variate from the *normal inverse gaussian*
    * distribution with parameters @f$\alpha@f$, @f$\beta@f$ = `beta`,
    * @f$\mu@f$ = `mu` and @f$\delta@f$, using generators `ig` and `ng`,
    * as described in eq. (
    * {@link REF_randvar_NormalInverseGaussianIGGen_nig2
    * nig2} ). The parameters @f$\alpha@f$ and @f$\delta@f$ are included
    * in generator `ig`.
    */
   public static double nextDouble (InverseGaussianGen ig, NormalGen ng,
                                    double beta, double mu) {
      return mynig (ig, ng, beta, mu);
   }
 

   public double nextDouble() {
      return mynig (genIG, genN, beta, mu);
   }


// >>>>>>>>>>>>>>>>>>>>  P R I V A T E     M E T H O D S   <<<<<<<<<<<<<<<<<<<

   private static double mynig (InverseGaussianGen ig, NormalGen ng,
                                double beta, double mu) {

      double y = ig.nextDouble ();
      double x = mu + beta*y + Math.sqrt(y)*ng.nextDouble ();
      return x;
   }


   protected void setParams (InverseGaussianGen ig, NormalGen ng,
                             double beta, double mu) {
      delta = Math.sqrt(ig.getLambda());
      gamma = delta / ig.getMu();
      alpha = Math.sqrt(gamma*gamma + beta*beta);
      setParams (alpha, beta, mu, delta);
      this.genN = ng;
      this.genIG = ig;
   }
}
