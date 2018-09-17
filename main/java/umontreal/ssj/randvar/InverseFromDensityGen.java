/*
 * Class:        InverseFromDensityGen
 * Description:  generator of random variates by numerical inversion of
                 an arbitrary continuous distribution when only the
                 probability density is known
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        June 2008
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
   import umontreal.ssj.functions.MathFunction;
   import umontreal.ssj.rng.RandomStream;
   import umontreal.ssj.probdist.ContinuousDistribution;
   import umontreal.ssj.probdist.InverseDistFromDensity;

/**
 * Implements a method for generating random variates by numerical inversion
 * of an *arbitrary continuous* distribution when only the probability
 * density is known @cite rDER10a&thinsp;. The cumulative probabilities (cdf)
 * are pre-computed by numerical quadrature of the density using
 * Gauss-Lobatto integration over suitably small intervals to satisfy the
 * required precision, and these values are kept in tables. Then the
 * algorithm uses polynomial interpolation over the tabulated values to get
 * the inverse cdf. The user can select the desired precision and the degree
 * of the interpolating polynomials.
 *
 * The algorithm may fail for some distributions for which the density
 * becomes infinite at a point (for ex. the Gamma and the Beta distributions
 * with @f$\alpha< 1@f$) if one requires too high a precision (a too small
 * `eps`, for ex. @f$\epsilon\sim10^{-15}@f$). However, it should work also
 * for continuous densities with finite discontinuities.
 *
 * While the setup time is relatively slow, the generation of random
 * variables is extremely fast and practically independent of the required
 * precision and of the specific distribution. The following table shows the
 * time needed (in seconds) to generate @f$10^8@f$ random numbers using
 * inversion from a given class, then the numerical inversion with
 * Gauss-Lobatto integration implemented here, and finally the speed ratios
 * between the two methods. The speed ratio is the speed of the latter over
 * the former. Thus for the beta distribution with parameters (5, 500),
 * generating random variables with the Gauss-Lobatto integration implemented
 * in this class is more than 1700 times faster than using inversion from the
 * `BetaDist` class. These tests were made on a machine with processor AMD
 * Athlon 4000, running Red Hat Linux, with clock speed at 2403 MHz.
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l bl br">Distribution</td>
 *   <td class="c bl br">Inversion</td>
 *   <td class="c bl br">Gauss-Lobatto</td>
 *   <td class="c bl br">speed ratio</td>
 * </tr><tr class="bt">
 *   <td class="l bl br">`NormalDist(10.5, 5)`</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;&ensp;9.19</td>
 *   <td class="c bl br">8.89</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;&ensp;&ensp;1.03</td>
 * </tr><tr>
 *   <td class="l bl br">`ExponentialDist(5)`</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;17.72</td>
 *   <td class="c bl br">8.82</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;&ensp;2.0</td>
 * </tr><tr>
 *   <td class="l bl br">`CauchyDist(10.5, 5)`</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;18.30</td>
 *   <td class="c bl br">8.81</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;&ensp;2.1</td>
 * </tr><tr>
 *   <td class="l bl br">`BetaSymmetricalDist(10.5)`</td>
 *   <td class="c bl br">&ensp;&ensp;242.80</td>
 *   <td class="c bl br">8.85</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;27.4</td>
 * </tr><tr>
 *   <td class="l bl br">`GammaDist(55)`</td>
 *   <td class="c bl br">&ensp;&ensp;899.50</td>
 *   <td class="c bl br">8.89</td>
 *   <td class="c bl br">&ensp;101</td>
 * </tr><tr>
 *   <td class="l bl br">`ChiSquareNoncentralDist(10.5, 5)`</td>
 *   <td class="c bl br">&ensp;5326.90</td>
 *   <td class="c bl br">8.85</td>
 *   <td class="c bl br">&ensp;602</td>
 * </tr><tr>
 *   <td class="l bl br">`BetaDist(5, 500)`</td>
 *   <td class="c bl br">15469.10</td>
 *   <td class="c bl br">8.86</td>
 *   <td class="c bl br">1746</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * The following table gives the time (in sec.) needed to create an object
 * (setup time) and to generate one random variable for this class compared
 * to the same for the inversion method specific to each class, and the
 * ratios of the times (init + one random variable) of the two methods. For
 * inversion, we initialized @f$10^8@f$ times; for this class, we initialized
 * @f$10^4@f$ times.
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l bl br">Distribution</td>
 *   <td class="c bl br">Inversion</td>
 *   <td class="c bl br">Gauss-Lobatto</td>
 *   <td class="c bl br">time ratio</td>
 * </tr><tr>
 *   <td class="l bl br"></td>
 *   <td class="c bl br">@f$10^8@f$ init</td>
 *   <td class="c bl br">@f$10^4@f$ init</td>
 *   <td class="c bl br">for 1 init</td>
 * </tr><tr class="bt">
 *   <td class="l bl br">`NormalDist(10.5, 5)`</td>
 *   <td class="c bl br">&ensp;&ensp;5.30</td>
 *   <td class="c bl br">&ensp;38.29</td>
 *   <td class="c bl br">26426</td>
 * </tr><tr>
 *   <td class="l bl br">`ExponentialDist(5)`</td>
 *   <td class="c bl br">&ensp;&ensp;3.98</td>
 *   <td class="c bl br">&ensp;27.05</td>
 *   <td class="c bl br">12466</td>
 * </tr><tr>
 *   <td class="l bl br">`CauchyDist(10.5, 5)`</td>
 *   <td class="c bl br">&ensp;&ensp;5.05</td>
 *   <td class="c bl br">&ensp;58.39</td>
 *   <td class="c bl br">25007</td>
 * </tr><tr>
 *   <td class="l bl br">`BetaSymmetricalDist(10.5)`</td>
 *   <td class="c bl br">&ensp;90.66</td>
 *   <td class="c bl br">&ensp;68.33</td>
 *   <td class="c bl br">&ensp;2049</td>
 * </tr><tr>
 *   <td class="l bl br">`GammaDist(55)`</td>
 *   <td class="c bl br">&ensp;13.15</td>
 *   <td class="c bl br">&ensp;58.34</td>
 *   <td class="c bl br">&ensp;&ensp;639</td>
 * </tr><tr>
 *   <td class="l bl br">`ChiSquareNoncentralDist(10.5, 5)`</td>
 *   <td class="c bl br">190.48</td>
 *   <td class="c bl br">248.98</td>
 *   <td class="c bl br">&ensp;&ensp;451</td>
 * </tr><tr>
 *   <td class="l bl br">`BetaDist(5, 500)`</td>
 *   <td class="c bl br">&ensp;63.60</td>
 *   <td class="c bl br">116.57</td>
 *   <td class="c bl br">&ensp;&ensp;&ensp;75</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * If only a few random variables are needed, then using this class is not
 * efficient because of the slow set-up. But if one wants to generate large
 * samples from the same distribution with fixed parameters, then this class
 * will be very efficient. The following table gives the number of random
 * variables generated beyond which, using this class will be worthwhile.
 *
 * <center>
 *
 * <table class="SSJ-table SSJ-has-hlines">
 * <tr class="bt">
 *   <td class="l bl br">Distribution</td>
 *   <td class="c bl br">number of generated variables</td>
 * </tr><tr class="bt">
 *   <td class="l bl br">`NormalDist(10.5, 5)`</td>
 *   <td class="c bl br">41665</td>
 * </tr><tr>
 *   <td class="l bl br">`ExponentialDist(5)`</td>
 *   <td class="c bl br">15266</td>
 * </tr><tr>
 *   <td class="l bl br">`CauchyDist(10.5, 5)`</td>
 *   <td class="c bl br">31907</td>
 * </tr><tr>
 *   <td class="l bl br">`BetaSymmetricalDist(10.5)`</td>
 *   <td class="c bl br">2814</td>
 * </tr><tr>
 *   <td class="l bl br">`GammaDist(55)`</td>
 *   <td class="c bl br">649</td>
 * </tr><tr>
 *   <td class="l bl br">`ChiSquareNoncentralDist(10.5, 5)`</td>
 *   <td class="c bl br">467</td>
 * </tr><tr>
 *   <td class="l bl br">`BetaDist(5, 500)`</td>
 *   <td class="c bl br">75</td>
 * </tr>
 * </table>
 *
 * </center>
 *
 * Thus, for example, if one needs to generate less than 15266 exponential
 * random variables, then using the `InverseFromDensityGen` class is not
 * wortwhile: it will be faster to use inversion from the `ExponentialGen`
 * class.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup randvar_general
 */
public class InverseFromDensityGen extends RandomVariateGen {

   /**
    * Creates a new generator for the *continuous* distribution `dis`,
    * using stream `stream`. `dis` must have a well-defined density
    * method; its other methods are unused. For a non-standard
    * distribution `dis`, the user may wish to set the left and the right
    * boundaries between which the density is non-zero by calling methods
    * umontreal.ssj.probdist.ContinuousDistribution.setXinf and
    * umontreal.ssj.probdist.ContinuousDistribution.setXsup of `dis`, for
    * better efficiency. Argument `xc` can be the mean, the mode or any
    * other @f$x@f$ for which the density is relatively large. The
    * @f$u@f$-resolution `eps` is the desired absolute error in the CDF,
    * and `order` is the degree of the Newton interpolating polynomial
    * over each interval. An `order` of 3 or 5, and an `eps` of
    * @f$10^{-6}@f$ to @f$10^{-12}@f$ are usually good choices.
    * Restrictions: @f$3 \le\mathtt{order} \le12@f$.
    */
   public InverseFromDensityGen (RandomStream stream,
                                 ContinuousDistribution dis,
                                 double xc, double eps, int order) {
      super (stream, null);
      dist = new InverseDistFromDensity (dis, xc, eps, order);
   }

   /**
    * Creates a new generator from the *continuous* probability density
    * `dens`. The left and the right boundaries of the density are `xleft`
    * and `xright` (the density is 0 outside the interval <tt>[xleft,
    * xright]</tt>). See the description of the other constructor.
    */
   public InverseFromDensityGen (RandomStream stream, MathFunction dens,
                                 double xc, double eps, int order,
                                 double xleft, double xright) {
      super (stream, null);
      dist = new InverseDistFromDensity (dens, xc, eps, order, xleft, xright);
   }

   /**
    * Generates a new random variate.
    */
   public double nextDouble() {
      return dist.inverseF (stream.nextDouble());
   }

   /**
    * Returns the `xc` given in the constructor.
    */
   public double getXc() {
      return ((InverseDistFromDensity)dist).getXc();
   }

   /**
    * Returns the @f$u@f$-resolution `eps`.
    */
   public double getEpsilon() {
      return ((InverseDistFromDensity)dist).getEpsilon();
   }

   /**
    * Returns the order of the interpolating polynomial.
    */
   public int getOrder() {
      return ((InverseDistFromDensity)dist).getOrder();
   }
}