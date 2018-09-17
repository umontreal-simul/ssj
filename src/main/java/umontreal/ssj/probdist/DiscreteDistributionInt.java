/*
 * Class:        DiscreteDistributionInt
 * Description:  discrete distributions over the integers
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
package umontreal.ssj.probdist;

/**
 * Classes implementing discrete distributions over the integers should inherit from this class. It
 * specifies the signatures of methods for computing the mass function (or probability) @f$p(x) =
 * P[X=x]@f$, distribution function @f$F(x)@f$, complementary distribution function
 * @f$\bar{F}(x)@f$, and inverse distribution function @f$F^{-1}(u)@f$, for a random
 *                   variable @f$X@f$ with a discrete distribution over the integers.
 *                   
 *                   *WARNING:* the complementary distribution function is defined as
 * @f$\bar{F}(j) = P[X \ge j]@f$ (for integers @f$j@f$, so that for discrete distributions in
 *               SSJ, @f$F(j) + \bar{F}(j) \neq1@f$ since both include the term @f$P[X = j]@f$.
 *
 *               The implementing classes provide both static and non-static methods to compute the
 *               above functions. The non-static methods require the creation of an object of class
 * @ref umontreal.ssj.probdist.DiscreteDistributionInt; all the non-negligible terms of the mass and
 *      distribution functions will be precomputed by the constructor and kept in arrays. Subsequent
 *      accesses will be very fast. The static methods do not require the construction of an object.
 *      These static methods are not specified in this abstract class because the number and types
 *      of their parameters depend on the distribution. When methods have to be called several times
 *      with the same parameters for the distributions, it is usually more efficient to create an
 *      object and use its non-static methods instead of the static ones. This trades memory for
 *      speed.
 *
 *      <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_general
 */

public abstract class DiscreteDistributionInt implements Distribution {

	/**
	 * Environment variable that determines what probability terms can be considered as negligible
	 * when building precomputed tables for distribution and mass functions. Probabilities smaller
	 * than `EPSILON` are not stored in the @ref umontreal.ssj.probdist.DiscreteDistribution objects
	 * (such as those of class @ref PoissonDist, etc.), but are computed directly each time they are
	 * needed (which should be very seldom). The default value is set to @f$10^{-16}@f$.
	 */
	public static double EPSILON = 1.0e-16;

	/*
	 * For better precision in the tails, we keep the cumulative probabilities (F) in cdf[x] for x
	 * <= xmed (i.e. cdf[x] is the sum off all the probabi- lities pdf[i] for i <= x), and the
	 * complementary cumulative probabilities (1 - F) in cdf[x] for x > xmed (i.e. cdf[x] is the sum
	 * off all the probabilities pdf[i] for i >= x).
	 */
	protected final static double EPS_EXTRA = 1.0e-6;
	protected double cdf[] = null;    // cumulative probabilities
	protected double pdf[] = null;    // probability terms or mass distribution
	protected int xmin = 0;           // pdf[x] < EPSILON for x < xmin
	protected int xmax = 0;           // pdf[x] < EPSILON for x > xmax

	// xmed is such that cdf[xmed] >= 0.5 and cdf[xmed - 1] < 0.5.
	protected int xmed = 0; /*
	                         * cdf[x] = F(x) for x <= xmed, and cdf[x] = bar_F(x) for x > xmed
	                         */
	protected int supportA = Integer.MIN_VALUE;
	protected int supportB = Integer.MAX_VALUE;

	/**
	 * Returns @f$p(x)@f$, the probability of @f$x@f$.
	 * 
	 * @param x
	 *            value at which the mass function must be evaluated
	 * @return the mass function evaluated at `x`
	 */
	public abstract double prob(int x);

	/**
	 * Returns the distribution function @f$F@f$ evaluated at @f$x@f$ (see (
	 * {@link REF_probdist_overview_eq_FDistDisc FDistDisc} )). Calls the {@link #cdf(int) cdf(int)}
	 * method.
	 * 
	 * @param x
	 *            value at which the distribution function must be evaluated
	 * @return the distribution function evaluated at `x`
	 */
	public double cdf(double x) {
		return cdf((int) x);
	}

	/**
	 * Returns the distribution function @f$F@f$ evaluated at @f$x@f$ (see (
	 * {@link REF_probdist_overview_eq_FDistDisc FDistDisc} )).
	 * 
	 * @param x
	 *            value at which the distribution function must be evaluated
	 * @return the distribution function evaluated at `x`
	 */
	public abstract double cdf(int x);

	/**
	 * Returns @f$\bar{F}(x)@f$, the complementary distribution function. Calls the
	 * {@link #barF(int) barF(int)} method.
	 * 
	 * @param x
	 *            value at which the complementary distribution function must be evaluated
	 * @return the complementary distribution function evaluated at `x`
	 */
	public double barF(double x) {
		return barF((int) x);
	}

	/**
	 * Returns @f$\bar{F}(x)@f$, the complementary distribution function. *See the WARNING above*.
	 * The default implementation returns `1.0 - cdf(x - 1)`, which is not accurate when @f$F(x)@f$
	 * is near 1.
	 * 
	 * @param x
	 *            value at which the complementary distribution function must be evaluated
	 * @return the complementary distribution function evaluated at `x`
	 */
	public double barF(int x) {
		return 1.0 - cdf(x - 1);
	}

	/**
	 * Returns the lower limit @f$x_a@f$ of the support of the probability mass function. The
	 * probability is 0 for all @f$x < x_a@f$.
	 * 
	 * @return @f$x@f$ lower limit of support
	 */
	public int getXinf() {
		return supportA;
	}

	/**
	 * Returns the upper limit @f$x_b@f$ of the support of the probability mass function. The
	 * probability is 0 for all @f$x > x_b@f$.
	 * 
	 * @return @f$x@f$ upper limit of support
	 */
	public int getXsup() {
		return supportB;
	}

	/**
	 * Returns the inverse distribution function @f$F^{-1}(u)@f$, where
	 * 
	 * @f$0\le u\le1@f$. Calls the `inverseFInt` method.
	 * @param u
	 *            value in the interval @f$(0,1)@f$ for which the inverse distribution function is
	 *            evaluated
	 * @return the inverse distribution function evaluated at `u`
	 *
	 * @exception IllegalArgumentException
	 *                if @f$u@f$ is not in the interval @f$(0,1)@f$
	 * @exception ArithmeticException
	 *                if the inverse cannot be computed, for example if it would give infinity in a
	 *                theoritical context
	 */
	public double inverseF(double u) {
		return inverseFInt(u);
	}

	/**
	 * Returns the inverse distribution function @f$F^{-1}(u)@f$, where
	 * 
	 * @f$0\le u\le1@f$. The default implementation uses binary search.
	 * @param u
	 *            value in the interval @f$(0,1)@f$ for which the inverse distribution function is
	 *            evaluated
	 * @return the inverse distribution function evaluated at `u`
	 *
	 * @exception IllegalArgumentException
	 *                if @f$u@f$ is not in the interval @f$(0,1)@f$
	 * @exception ArithmeticException
	 *                if the inverse cannot be computed, for example if it would give infinity in a
	 *                theoritical context
	 */
	public int inverseFInt(double u) {
		int i, j, k;

		if (u < 0.0 || u > 1.0)
			throw new IllegalArgumentException("u is not in [0,1]");
		if (u <= 0.0)
			return supportA;
		if (u >= 1.0)
			return supportB;

		// Remember: the upper part of cdf contains the complementary distribu-
		// tion for xmed < s <= xmax, and the lower part of cdf the
		// distribution for xmin <= x <= xmed

		if (u <= cdf[xmed - xmin]) {
			// In the lower part of cdf
			if (u <= cdf[0])
				return xmin;
			i = 0;
			j = xmed - xmin;
			while (i < j) {
				k = (i + j) / 2;
				if (u > cdf[k])
					i = k + 1;
				else
					j = k;
			}
		} else {
			// In the upper part of cdf
			u = 1 - u;
			if (u < cdf[xmax - xmin])
				return xmax;

			i = xmed - xmin + 1;
			j = xmax - xmin;
			while (i < j) {
				k = (i + j) / 2;
				if (u < cdf[k])
					i = k + 1;
				else
					j = k;
			}
			i--;
		}

		return i + xmin;
	}

}