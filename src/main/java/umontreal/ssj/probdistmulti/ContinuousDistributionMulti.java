/*
 * Class:        ContinuousDistributionMulti
 * Description:  mother class for continuous multidimensional distributions 
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

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