/*
 * Class:        GeometricNormalInverseGaussianProcess
 * Description:  
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
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.randvar.*;

/**
 * The geometric normal inverse gaussian (GNIG) process is the exponentiation
 * of a  @ref NormalInverseGaussianProcess :
 * @f[
 *   S(t) = S_0 \exp\left[ (r-\omega_{RN})t + \mbox{NIG}(t;\alpha,\beta,\mu,\delta) \right],
 * @f]
 * where @f$r@f$ is the interest rate. It is a strictly positive process,
 * which is useful in finance. There is also a neutral correction in the
 * exponential, @f$\omega_{RN}= \mu+
 * \delta\gamma-\delta\sqrt{\alpha^2-(1+\beta)^2}@f$, which takes into
 * account the market price of risk. The underlying NIG process must start at
 * zero, NIG@f$(t_0) = 0 @f$ and the initial time should also be set to zero,
 * @f$t_0 = 0@f$, both for the NIG and GNIG.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class GeometricNormalInverseGaussianProcess extends
                                                   GeometricLevyProcess {

   /**
    * Constructs a new <tt>GeometricNormalInverseGaussianProcess</tt>. The
    * parameters of the NIG process will be overwritten by the parameters
    * given to the GNIG, with the initial value of the NIG set to 0. The
    * observation times of the NIG will also be changed to those of the
    * GNIG.
    */
   public GeometricNormalInverseGaussianProcess (
                                        double s0, double muGeom,
                                        double alpha, double beta,
                                        double mu, double delta,
                                        RandomStream streamBrownian,
                                        NormalInverseGaussianProcess nigP) {
        levyProcess = nigP;
        ((NormalInverseGaussianProcess)levyProcess).setParams(0.0, alpha, beta, mu, delta);
        this.x0 = s0;
        this.muGeom = muGeom;
        omegaRiskNeutralCorrection =
            mu + delta*((NormalInverseGaussianProcess)levyProcess).getGamma() -
            delta*Math.sqrt(alpha*alpha - (1.+beta)*(1.+beta));
    }

   /**
    * Constructs a new `GeometricNormalInverseGaussianProcess`. The
    * process `igP` will be used internally by the underlying
    * @ref NormalInverseGaussianProcess.
    */
   public GeometricNormalInverseGaussianProcess (
                                        double s0, double muGeom,
                                        double alpha, double beta,
                                        double mu, double delta,
                                        RandomStream streamBrownian,
                                        InverseGaussianProcess igP) {

        levyProcess = new NormalInverseGaussianProcess (0.0, alpha, beta, mu,
                                                        delta, streamBrownian, igP);
        this.x0 = s0;
        this.muGeom = muGeom;
        omegaRiskNeutralCorrection =
            mu + delta*((NormalInverseGaussianProcess)levyProcess).getGamma() -
            delta*Math.sqrt(alpha*alpha - (1.+beta)*(1.+beta));
    }

   /**
    * Constructs a new `GeometricNormalInverseGaussianProcess`. The drift
    * of the geometric term, `muGeom`, is usually the interest rate
    * @f$r@f$. `s0` is the initial value of the process and the other four
    * parameters are the parameters of the underlying
    * @ref NormalInverseGaussianProcess process.
    */
   public GeometricNormalInverseGaussianProcess (
                                        double s0, double muGeom,
                                        double alpha, double beta,
                                        double mu, double delta,
                                        RandomStream streamBrownian,
                                        RandomStream streamNIG1,
                                        RandomStream streamNIG2,
                                        String igType) {
        levyProcess = new NormalInverseGaussianProcess (0.0, alpha, beta, mu, delta,
                                                        streamBrownian, streamNIG1,
                                                        streamNIG2, igType);
        this.x0 = s0;
        this.muGeom = muGeom;
        omegaRiskNeutralCorrection =
            mu + delta*((NormalInverseGaussianProcess)levyProcess).getGamma() -
            delta*Math.sqrt(alpha*alpha - (1.+beta)*(1.+beta));
    }

   /**
    * Constructs a new `GeometricNormalInverseGaussianProcess`. The String
    * `igType` corresponds to the type of  @ref InverseGaussianProcess
    * that will be used by the underlying
    * @ref NormalInverseGaussianProcess. All
    * @ref umontreal.ssj.rng.RandomStream ’s used to generate the
    * underlying  @ref NormalInverseGaussianProcess and its underlying
    * @ref InverseGaussianProcess are set to the same given `streamAll`.
    */
   public GeometricNormalInverseGaussianProcess (
                                        double s0, double muGeom,
                                        double alpha, double beta,
                                        double mu, double delta,
                                        RandomStream streamAll,
                                        String igType) {
        this(s0,muGeom,alpha,beta,mu,delta,streamAll,streamAll,streamAll,igType);
    }

}
