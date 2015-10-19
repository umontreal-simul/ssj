/*
 * Class:        VarianceGammaProcessDiffPCABridge
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @authors      Jean-Sébastien Parent & Maxime Dion
 * @since        2008

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
package umontreal.ssj.stochprocess;
import umontreal.ssj.rng.*;

/**
 * Same as  @ref VarianceGammaProcessDiff, but the two inner
 * @ref GammaProcess ’es are of the type PCABridge. Also,
 * `generatePath(double[] uniform01)` distributes the lowest coordinates
 * uniforms to the inner  @ref GammaProcessPCABridge according to their
 * eigenvalues.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class VarianceGammaProcessDiffPCABridge extends
                                               VarianceGammaProcessDiffPCA {

   /**
    * Constructs a new  @ref VarianceGammaProcessDiffPCABridge with
    * parameters @f$\theta= \mathtt{theta}@f$, @f$\sigma=
    * \mathtt{sigma}@f$, @f$\nu= \mathtt{nu}@f$ and initial value
    * @f$S(t_0) = \mathtt{s0}@f$. There is only one
    * @ref umontreal.ssj.rng.RandomStream here which is used for the two
    * inner  @ref GammaProcessPCABridge ’s. The other parameters are set
    * as in  @ref VarianceGammaProcessDiff.
    */
   public VarianceGammaProcessDiffPCABridge (double s0, double theta,
                                             double sigma, double nu,
                                             RandomStream stream) {
     super(s0, theta, sigma, nu,
	  new GammaProcessPCABridge (0.0, 1.0, 1.0, stream),
	  new GammaProcessPCABridge (0.0, 1.0, 1.0, stream));
    // Params mu, nu of the 2 gamma processes are redefined in init()
    // which will be called after a call to 'setObservTimes'
}

}