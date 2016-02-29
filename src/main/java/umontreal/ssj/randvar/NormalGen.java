/*
 * Class:        NormalGen
 * Description:  random variates generator from the normal distribution
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
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;

/**
 * This class implements methods for generating random variates from the
 * *normal* distribution @f$N(\mu, \sigma)@f$. It has mean @f$\mu@f$ and
 * variance @f$\sigma^2@f$, where @f$\sigma>0@f$. Its density function is
 * @anchor REF_randvar_NormalGen_eq_fnormal
 * @f[
 *   f(x) = \frac{1}{\sqrt{2\pi}\sigma} e^{(x-\mu)^2/(2\sigma^2)} \tag{fnormal}
 * @f]
 * The `nextDouble` method simply calls `inverseF` on the distribution.
 *
 * The following table gives the CPU time needed to generate @f$10^8@f$
 * standard normal random variates using the different implementations
 * available in SSJ. The first time is for a generator object (non-static
 * method), and the second time is for the static method where no object is
 * created. These tests were made on a machine with processor AMD Athlon
 * 4000, running Red Hat Linux, with clock speed at 2403 MHz. The static
 * method `nextDouble()` for `NormalBoxMullerGen` and `NormalPolarGen` uses
 * only one number out of two that are generated; thus they are twice slower
 * than the non-static method.
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l bl br">Generator</td>
 *   <td class="c bl br">time in seconds</td>
 *   <td class="c bl br">time in seconds</td>
 * </tr><tr>
 *   <td class="l bl br"></td>
 *   <td class="c bl br">(object)</td>
 *   <td class="c bl br">(static)</td>
 * </tr><tr class="bt">
 *   <td class="l bl br">`NormalGen`</td>
 *   <td class="c bl br">7.67</td>
 *   <td class="c bl br">7.72</td>
 * </tr><tr>
 *   <td class="l bl br">`NormalACRGen`</td>
 *   <td class="c bl br">4.71</td>
 *   <td class="c bl br">4.76</td>
 * </tr><tr>
 *   <td class="l bl br">`NormalBoxMullerGen`</td>
 *   <td class="c bl br">16.07</td>
 *   <td class="c bl br">31.45</td>
 * </tr><tr>
 *   <td class="l bl br">`NormalPolarGen`</td>
 *   <td class="c bl br">7.31</td>
 *   <td class="c bl br">13.74</td>
 * </tr><tr>
 *   <td class="l bl br">`NormalKindermannRamageGen`</td>
 *   <td class="c bl br">5.38</td>
 *   <td class="c bl br">5.34</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_continuous
 */
public class NormalGen extends RandomVariateGen {
   protected double mu;
   protected double sigma = -1.0;

   /**
    * Creates a normal random variate generator with mean `mu` and
    * standard deviation `sigma`, using stream `s`.
    */
   public NormalGen (RandomStream s, double mu, double sigma) {
      super (s, new NormalDist(mu, sigma));
      setParams (mu, sigma);
   }

   /**
    * Creates a standard normal random variate generator with mean `0` and
    * standard deviation `1`, using stream `s`.
    */
   public NormalGen (RandomStream s) {
      this (s, 0.0, 1.0);
   }

   /**
    * Creates a random variate generator for the normal distribution
    * `dist` and stream `s`.
    */
   public NormalGen (RandomStream s, NormalDist dist) {
      super (s, dist);
      if (dist != null)
         setParams (dist.getMu(), dist.getSigma());
   }

   /**
    * Generates a variate from the normal distribution with parameters
    * @f$\mu= @f$&nbsp;`mu` and @f$\sigma= @f$&nbsp;`sigma`, using
    * stream `s`.
    */
   public static double nextDouble (RandomStream s, double mu, double sigma) {
      return NormalDist.inverseF (mu, sigma, s.nextDouble());
   }

   /**
    * Returns the parameter @f$\mu@f$ of this object.
    */
   public double getMu() {
      return mu;
   }

   /**
    * Returns the parameter @f$\sigma@f$ of this object.
    */
   public double getSigma() {
      return sigma;
   }

   /**
    * Sets the parameters @f$\mu@f$ and @f$\sigma@f$ of this object.
    */
   protected void setParams (double mu, double sigma) {
      if (sigma <= 0)
         throw new IllegalArgumentException ("sigma <= 0");
      this.mu = mu;
      this.sigma = sigma;
   }
}