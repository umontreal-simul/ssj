/*
 * Class:        KernelDensityVarCorrectGen
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
 * This class is a variant of  @ref KernelDensityGen, but with a rescaling of
 * the empirical distribution so that the variance of the density used to
 * generate the random variates is equal to the empirical variance, as
 * suggested by @cite tSIL86a&thinsp;.
 *
 * Let @f$\bar{x}_n@f$ and @f$s_n^2@f$ be the sample mean and sample variance
 * of the observations. The distance between each generated random variate
 * and the sample mean @f$\bar{x}_n@f$ is multiplied by the correcting factor
 * @f$1/\sigma_e@f$, where @f$\sigma_e^2 = 1 + (h\sigma_k/s_n)^2@f$. The
 * constant @f$\sigma_k^2@f$ must be passed to the constructor. Its value
 * can be found in Table&nbsp;
 * {@link REF_randvar_KernelDensityGen_tab_kernels kernels}
 * for some popular kernels.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class KernelDensityVarCorrectGen extends KernelDensityGen {

   protected double sigmak2;   // Value of sigma_k^2.
   protected double mean;      // Sample mean of the observations.
   protected double invSigmae; // 1 / sigma_e.

   /**
    * Creates a new generator for a kernel density estimated from the
    * observations given by the empirical distribution `dist`, using
    * stream `s` to select the observations, generator `kGen` to generate
    * the added noise from the kernel density, bandwidth `h`, and
    * @f$\sigma_k^2 =@f$ `sigmak2` used for the variance correction.
    */
   public KernelDensityVarCorrectGen (RandomStream s, EmpiricalDist dist,
                                      RandomVariateGen kGen, double h,
                                      double sigmak2) {
      super (s, dist, kGen, h);
      this.sigmak2 = sigmak2;
      mean = dist.getSampleMean();
      double var = dist.getSampleVariance();
      invSigmae = 1.0 / Math.sqrt (1.0 + h * h * sigmak2 / var);
   }

   /**
    * This constructor uses a gaussian kernel and the default bandwidth
    * suggested in Table&nbsp;
    * {@link REF_randvar_KernelDensityGen_tab_kernels
    * kernels} for the gaussian distribution.
    */
   public KernelDensityVarCorrectGen (RandomStream s, EmpiricalDist dist,
                                      NormalGen kGen) {
      this (s, dist, kGen, 0.77639 * getBaseBandwidth (dist), 1.0);
   }


   public void setBandwidth (double h) {
      if (h < 0)
         throw new IllegalArgumentException ("h < 0");
      bandwidth = h;
      double var = ((EmpiricalDist) dist).getSampleVariance();
      invSigmae = 1.0 / Math.sqrt (1.0 + h * h * sigmak2 / var);
   }

   public double nextDouble() {
      double x = mean + invSigmae * (dist.inverseF (stream.nextDouble())
                  - mean + bandwidth * kernelGen.nextDouble());
      if (positive)
         return Math.abs (x);
      else
         return x;
   }
}