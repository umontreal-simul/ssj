/*
 * Class:        ContinuousDistributionMulti
 * Description:  mother class for continuous multidimensional distributions 
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
package umontreal.ssj.probdistmulti;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.util.Num;

/**
 * Classes implementing continuous multi-dimensional distributions should
 * inherit from this class. Such distributions are characterized by a
 * *density* function @f$f(x_1, x_2, …, x_d)@f$; thus the signature of a
 * `density` method is supplied here. All array indices start at 0.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_general
 */
public abstract class ContinuousDistributionMulti {
   protected int dimension;

/**
 * Returns @f$f(x_1, x_2, …, x_d)@f$, the probability density of @f$X@f$
 * evaluated at the point @f$x@f$, where @f$x = \{x_1, x_2, …, x_d\}@f$. The
 * convention is that @f$\mathtt{x[i-1]} = x_i@f$.
 *  @param x            value at which the density is evaluated
 *  @return density function evaluated at `x`
 */
public abstract double density (double[] x);

   /**
    * Returns the dimension @f$d@f$ of the distribution.
    */
   public int getDimension() {
      return dimension;
   }

   /**
    * Returns the mean vector of the distribution, defined as @f$\mu_i =
    * E[X_i]@f$.
    */
   public abstract double[] getMean();

   /**
    * Returns the variance-covariance matrix of the distribution, defined
    * as<br>@f$\sigma_{ij} = E[(X_i - \mu_i)(X_j - \mu_j)]@f$.
    */
   public abstract double[][] getCovariance();

   /**
    * Returns the correlation matrix of the distribution, defined as
    * @f$\rho_{ij} = \sigma_{ij}/\sqrt{\sigma_{ii}\sigma_{jj}}@f$.
    */
   public abstract double[][] getCorrelation();

}