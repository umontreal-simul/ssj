/*
 * Class:        VarianceGammaProcessDiffPCASymmetricalBridge
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @authors       Jean-Sébastien Parent-Chartier and Maxime Dion 
 * @since        2008
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

/**
 * Same as  @ref VarianceGammaProcessDiff, but the two inner
 * @ref GammaProcess ’es are of the PCASymmetricalBridge type. Also,
 * `generatePath(double[] uniform01)` distributes the lowest coordinates
 * uniforms to the inner GammaProcessPCA according to their eigenvalues.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class VarianceGammaProcessDiffPCASymmetricalBridge extends
             VarianceGammaProcessDiffPCA {

/**
 * Constructs a new  @ref VarianceGammaProcessDiffPCASymmetricalBridge with
 * parameters @f$\theta= \mathtt{theta}@f$, @f$\sigma= \mathtt{sigma}@f$,
 * @f$\nu= \mathtt{nu}@f$ and initial value @f$S(t_0) = \mathtt{s0}@f$.
 * There is only one  @ref umontreal.ssj.rng.RandomStream here which is used
 * for the two inner  @ref GammaProcessPCASymmetricalBridge ’s. The other
 * parameters are set as in  @ref VarianceGammaProcessDiff.
 */
public VarianceGammaProcessDiffPCASymmetricalBridge (
                                               double s0, double theta,
                                               double sigma, double nu,
                                               RandomStream stream) {
    super(s0, theta, sigma, nu, new GammaProcessPCASymmetricalBridge (0.0, 1.0, 1.0, stream),
	  new GammaProcessPCASymmetricalBridge (0.0, 1.0, 1.0, stream));
    // Params mu, nu of the 2 gamma processes are redefined in init()
    // which will be called after a call to 'setObservTimes'
}

}