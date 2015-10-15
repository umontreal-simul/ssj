/*
 * Class:        DiscreteDistributionIntMulti
 * Description:  mother class for discrete distributions over the integers
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

/**
 * Classes implementing multi-dimensional discrete distributions over the
 * integers should inherit from this class. It specifies the signature of
 * methods for computing the mass function (or probability) @f$p(x_1, x_2, …,
 * x_d) = P[X_1 = x_1, X_2 = x_2, …, X_d = x_d]@f$ and the cumulative
 * probabilities for a random vector @f$X@f$ with a discrete distribution
 * over the integers.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdistmulti_general
 */
public abstract class DiscreteDistributionIntMulti {
   protected int dimension;

/**
 * Returns the probability mass function @f$p(x_1, x_2, …, x_d)@f$, which
 * should be a real number in @f$[0,1]@f$.
 *  @param x            value at which the mass function must be evaluated
 *  @return the mass function evaluated at `x`
 */
public abstract double prob (int[] x);

   /**
    * Computes the cumulative probability function @f$F@f$ of the
    * distribution evaluated at `x`, assuming the lowest values start at
    * 0, i.e. computes
    * @f[
    *   F (x_1, x_2, …, x_d) = \sum_{s_1=0}^{x_1} \sum_{s_2=0}^{x_2} \cdots\sum_{s_d=0}^{x_d} p(s_1, s_2, …, s_d).
    * @f]
    * Uses the naive implementation, is very inefficient and may
    * underflows.
    */
   public double cdf (int x[]) {
      int is[] = new int[x.length];
      for (int i = 0; i < is.length; i++)
         is[i] = 0;

      boolean end = false;
      double sum = 0.0;
      int j;
      while (!end) {
         sum += prob (is);
         is[0]++;

         if (is[0] > x[0]) {
            is[0] = 0;
            j = 1;
            while (j < x.length && is[j] == x[j])
               is[j++] = 0;

            if (j == x.length)
               end = true;
            else
               is[j]++;
         }
      }

      return sum;
   }

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