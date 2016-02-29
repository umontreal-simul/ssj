/*
 * Class:        KernelDensityGen
 * Description:  random variate generators for distributions obtained via
                 kernel density estimation methods
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
package umontreal.ssj.randvar;
import umontreal.ssj.probdist.*;
import umontreal.ssj.rng.RandomStream;

/**
 * This class implements random variate generators for distributions obtained
 * via *kernel density estimation* methods from a set of @f$n@f$ individual
 * observations @f$x_1,…,x_n@f$ @cite tDEV85a, @cite rDEV86a, @cite rHOR00a,
 * @cite rHOR04a, @cite tSIL86a&thinsp;. The basic idea is to center a copy
 * of the same symmetric density at each observation and take an equally
 * weighted mixture of the @f$n@f$ copies as an estimator of the density from
 * which the observations come. The resulting kernel density has the general
 * form
 * @f[
 *   f_n(x) = \frac{1}{nh} \sum_{i=1}^n k((x-x_i)/h),
 * @f]
 * where @f$k@f$ is a fixed pre-selected density called the *kernel* and
 * @f$h@f$ is a positive constant called the *bandwidth* or *smoothing
 * factor*. A difficult practical issue is the selection of @f$k@f$ and
 * @f$h@f$. Several approaches have been proposed for that; see, e.g.,
 * @cite tBER94a, @cite tCAO94a, @cite rHOR04a, @cite tSIL86a&thinsp;.
 *
 * The constructor of a generator from a kernel density requires a random
 * stream @f$s@f$, the @f$n@f$ observations in the form of an empirical
 * distribution, a random variate generator for the kernel density @f$k@f$,
 * and the value of the bandwidth @f$h@f$. The random variates are then
 * generated as follows: select an observation @f$x_I@f$ at random, by
 * inversion, using stream @f$s@f$, then generate random variate @f$Y@f$ with
 * the generator provided for the density @f$k@f$, and return @f$x_I + hY@f$.
 *
 * A simple formula for the bandwidth, suggested in @cite tSIL86a,
 * @cite rHOR04a&thinsp;, is @f$h = \alpha_k h_0@f$, where
 * @anchor REF_randvar_KernelDensityGen_eq_bandwidth0
 * @f[
 *   h_0 = 1.36374 \min(s_n, q / 1.34) n^{-1/5}, \tag{bandwidth0}
 * @f]
 * @f$s_n@f$ and @f$q@f$ are the empirical standard deviation and the
 * interquartile range of the @f$n@f$ observations, and @f$\alpha_k@f$ is a
 * constant that depends on the type of kernel @f$k@f$. It is defined by
 * @f[
 *   \alpha_k = \left(\sigma_k^{-4} \int_{-\infty}^{\infty}k(x)dx \right)^{1/5}
 * @f]
 * where @f$\sigma_k@f$ is the standard deviation of the density @f$k@f$.
 * The static method  #getBaseBandwidth permits one to compute @f$h_0@f$ for
 * a given empirical distribution.
 *
 * <center>
 * <center>Some suggested kernels</center>
 * @anchor REF_randvar_KernelDensityGen_tab_kernels
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l">name</td>
 *   <td class="l">constructor</td>
 *   <td class="c">@f$\alpha_k@f$</td>
 *   <td class="c">@f$\sigma_k^2@f$</td>
 *   <td class="c">efficiency</td>
 * </tr><tr class="bt">
 *   <td class="l">Epanechnikov</td>
 *   <td class="l">`BetaSymmetricalDist(2, -1, 1)`</td>
 *   <td class="c">1.7188</td>
 *   <td class="c">1/5</td>
 *   <td class="c">1.000</td>
 * </tr><tr>
 *   <td class="l">triangular</td>
 *   <td class="l">`TriangularDist(-1, 1, 0)`</td>
 *   <td class="c">1.8882</td>
 *   <td class="c">1/6</td>
 *   <td class="c">0.986</td>
 * </tr><tr>
 *   <td class="l">Gaussian</td>
 *   <td class="l">`NormalDist()`</td>
 *   <td class="c">0.7764</td>
 *   <td class="c">1</td>
 *   <td class="c">0.951</td>
 * </tr><tr>
 *   <td class="l">boxcar</td>
 *   <td class="l">`UniformDist(-1, 1)`</td>
 *   <td class="c">1.3510</td>
 *   <td class="c">1/3</td>
 *   <td class="c">0.930</td>
 * </tr><tr>
 *   <td class="l">logistic</td>
 *   <td class="l">`LogisticDist()`</td>
 *   <td class="c">0.4340</td>
 *   <td class="c">3.2899</td>
 *   <td class="c">0.888</td>
 * </tr><tr>
 *   <td class="l">Student-t(3)</td>
 *   <td class="l">`StudentDist(3)`</td>
 *   <td class="c">0.4802</td>
 *   <td class="c">3</td>
 *   <td class="c">0.674</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * Table&nbsp; {@link REF_randvar_KernelDensityGen_tab_kernels
 * kernels} gives the precomputed values of @f$\sigma_k@f$ and
 * @f$\alpha_k@f$ for selected (popular) kernels. The values are taken from
 * @cite rHOR04a&thinsp;. The second column gives the name of a function (in
 * this package) that constructs the corresponding distribution. The
 * *efficiency* of a kernel is defined as the ratio of its mean integrated
 * square error over that of the Epanechnikov kernel, which has optimal
 * efficiency and corresponds to the beta distribution with parameters
 * @f$(2,2)@f$ over the interval @f$(-1,1)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class KernelDensityGen extends RandomVariateGen {

   protected RandomVariateGen kernelGen;
   protected double bandwidth;
   protected boolean positive;   // If we want positive reflection.

   /**
    * Creates a new generator for a kernel density estimated from the
    * observations given by the empirical distribution `dist`, using
    * stream `s` to select the observations, generator `kGen` to generate
    * the added noise from the kernel density, and bandwidth `h`.
    */
   public KernelDensityGen (RandomStream s, EmpiricalDist dist,
                            RandomVariateGen kGen, double h) {
      super (s, dist);
      if (h < 0.0)
         throw new IllegalArgumentException ("h < 0");
      if (kGen == null)
         throw new IllegalArgumentException ("kGen == null");
      kernelGen = kGen;
      bandwidth = h;
   }

   /**
    * This constructor uses a gaussian kernel and the default bandwidth
    * @f$h = \alpha_k h_0@f$ with the @f$\alpha_k@f$ suggested in
    * Table&nbsp;
    * {@link REF_randvar_KernelDensityGen_tab_kernels
    * kernels} for the gaussian distribution. This kernel has an
    * efficiency of 0.951.
    */
   public KernelDensityGen (RandomStream s, EmpiricalDist dist,
                            NormalGen kGen) {
      this (s, dist, kGen, 0.77639 * getBaseBandwidth (dist));
   }

   /**
    * @name Kernel selection and parameters
    * @{
    */

   /**
    * Computes and returns the value of @f$h_0@f$ in (
    * {@link REF_randvar_KernelDensityGen_eq_bandwidth0
    * bandwidth0} ).
    */
   public static double getBaseBandwidth (EmpiricalDist dist) {
      double r = dist.getInterQuartileRange() / 1.34;
      double sigma = dist.getSampleStandardDeviation();
      if (sigma < r) r = sigma;
      return (1.36374 * r / Math.exp (0.2 * Math.log (dist.getN())));
   }

   /**
    * Sets the bandwidth to `h`.
    */
   public void setBandwidth (double h) {
      if (h < 0)
         throw new IllegalArgumentException ("h < 0");
      bandwidth = h;
   }

   /**
    * After this method is called with `true`, the generator will produce
    * only positive values, by using the *reflection method*: replace all
    * negative values by their *absolute values*. That is,  #nextDouble
    * will return @f$|x|@f$ if @f$x@f$ is the generated variate. The
    * mecanism is disabled when the method is called with `false`.
    */
   public void setPositiveReflection (boolean reflect) {
      positive = reflect;
   }


   public double nextDouble() {
      double x = (dist.inverseF (stream.nextDouble())
                  + bandwidth * kernelGen.nextDouble());
      if (positive)
         return Math.abs (x);
      else
         return x;
   }
}

/**
 * @}
 */