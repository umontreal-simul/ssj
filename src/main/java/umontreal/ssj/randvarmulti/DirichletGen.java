/*
 * Class:        DirichletGen
 * Description:  multivariate generator for a Dirichlet distribution
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
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
package umontreal.ssj.randvarmulti;

import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.randvar.RandomVariateGen;
import umontreal.ssj.randvar.GammaAcceptanceRejectionGen;
import umontreal.ssj.rng.RandomStream;

/**
 * Extends  @ref RandomMultivariateGen for a *Dirichlet*
 * @cite tJOH72a&thinsp; distribution. This distribution uses the parameters
 * @f$\alpha_1,…,\alpha_k@f$, and has density
 * @f[
 *   f(x_1,…,x_k) = \frac{\Gamma(\alpha_0)\prod_{i=1}^k x_i^{\alpha_i - 1}}{\prod_{i=1}^k \Gamma(\alpha_i)}
 * @f]
 * where @f$\alpha_0=\sum_{i=1}^k\alpha_i@f$.
 *
 * Here, the successive coordinates of the Dirichlet vector are generated
 * @remark **Pierre:** How?
 *
 *  via the class  @ref umontreal.ssj.randvar.GammaAcceptanceRejectionGen in
 * package `randvar`, using the same stream for all the uniforms.
 *
 * Note: when the shape parameters @f$\alpha_i@f$ are all very small, the
 * results may lose some numerical precision. For example, the value of the
 * density function of the Dirichlet multivariate may return 0. Also, the
 * generated @f$x_i@f$ will often have one variate equals to 1, and all
 * others set at (or near) 0.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class DirichletGen extends RandomMultivariateGen {
   private GammaAcceptanceRejectionGen[] ggens;
   private double[] alphas;

   // use log gamma if all alphas are smaller than this threshold
   private static final double ALPHA_THRESHOLD = 0.1; 

   // determine if the Dirichlet multivariate must be generated using the
   // log gamma, when the shape parameter alphas are all small.
   private boolean useLogGamma = false;

   /**
    * Constructs a new Dirichlet generator with parameters
    * @f$\alpha_{i+1}=@f$&nbsp;`alphas[i]`, for @f$i=0,…,k-1@f$, and the
    * stream `stream`.
    *  @param stream       the random number stream used to generate
    *                      uniforms.
    *  @param alphas       the @f$\alpha_i@f$ parameters of the generated
    *                      distribution.
    *  @exception IllegalArgumentException if one @f$\alpha_k@f$ is
    * negative or 0.
    *  @exception NullPointerException if any argument is `null`.
    */
   public DirichletGen (RandomStream stream, double[] alphas) {
      if (stream == null)
         throw new NullPointerException ("stream is null");
      this.stream = stream;
      dimension = alphas.length;
      ggens = new GammaAcceptanceRejectionGen[alphas.length];
      for (int k = 0; k < alphas.length; k++)
         ggens[k] = new GammaAcceptanceRejectionGen
            (stream, new GammaDist (alphas[k], 1.0/2.0));

      this.alphas = new double[alphas.length];
      System.arraycopy(alphas, 0, this.alphas, 0, alphas.length);
      useLogGamma = canUseLogGamma(alphas);
   }

   /**
    * Returns the @f$\alpha_{i+1}@f$ parameter for this Dirichlet
    * generator.
    *  @param i            the index of the parameter.
    *  @return the value of the parameter.
    *
    *  @exception ArrayIndexOutOfBoundsException if `i` is negative or
    * greater than or equal to  #getDimension.
    */
   public double getAlpha (int i) {
      return alphas[i];
   }

   /**
    * Generates a new point from the Dirichlet distribution with
    * parameters `alphas`, using the stream `stream`. The generated values
    * are placed into `p`.
    *  @param stream       the random number stream used to generate the
    *                      uniforms.
    *  @param alphas       the @f$\alpha_i@f$ parameters of the
    *                      distribution, for @f$i=1,…,k@f$.
    *  @param p            the array to be filled with the generated
    *                      point.
    */
   public static void nextPoint (RandomStream stream, double[] alphas,
                                 double[] p) {
      if (canUseLogGamma(alphas)) {
         nextPointUsingLog(stream, alphas, p);
         return;
      }

      double total = 0;
      for (int i = 0; i < alphas.length; i++) {
         p[i] = GammaAcceptanceRejectionGen.nextDouble
            (stream, stream, alphas[i], 1.0/2.0);
         total += p[i];
      }
      for (int i = 0; i < alphas.length; i++)
         p[i] /= total;
   }

   /**
    * Generates a point from the Dirichlet distribution.
    *  @param p            the array to be filled with the generated
    *                      point.
    */
   public void nextPoint (double[] p) {
      if (useLogGamma) {
         nextPointUsingLog(stream, alphas, p);
         return;
      }

      int n = ggens.length;
      double total = 0;
      for (int i = 0; i < n; i++) {
         p[i] = ggens[i].nextDouble();
         total += p[i];
      }
      for (int i = 0; i < n; i++)
         p[i] /= total;
   }

   /**
    * Generates Dirichlet multivariate using the log of the gammas.
    * This should be used if the shape parameters are very small.
    */
   private static void nextPointUsingLog (RandomStream stream, double[] alphas,
                                 double[] p) {
      double total = 0;
      double[] log = new double[p.length]; // contains the log value

      // get the log value of the gammas
      for (int i = 0; i < alphas.length; i++) {
         log[i] = GammaAcceptanceRejectionGen.nextDoubleLog
            (stream, stream, alphas[i], 1.0/2.0);
      }

      // find the ratio for each p
      for (int i = 0; i < alphas.length; i++) {
          total = 0;
          for (int j = 0; j < alphas.length; j++)
             total += Math.exp(log[j]-log[i]);
          p[i] = 1.0 / total;
      }
   }

   /** 
    * Checks if all alphas are smaller than the threshold.
    * If it returns true, then the multivariate can be generated
    * using the log gammas.
    */
   private static boolean canUseLogGamma(double[] alphas) {
      for (int i = 0; i < alphas.length; i++)
         if (alphas[i] >= ALPHA_THRESHOLD)
            return false;
      return true;
   }

}