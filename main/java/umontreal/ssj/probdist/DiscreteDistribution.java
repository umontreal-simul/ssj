/*
 * Class:        DiscreteDistribution
 * Description:  discrete distributions over a set of real numbers
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

import java.util.Formatter;
import java.util.Locale;

/**
 * This class implements discrete distributions over a *finite set of real numbers* (also over
 * *integers* as a particular case). We assume that the random variable @f$X@f$ of interest can take
 * one of the @f$n@f$ values
 * 
 * @f$x_0 < \cdots< x_{n-1}@f$, which *must be sorted* by increasing order.
 * @f$X@f$ can take the value @f$x_k@f$ with probability @f$p_k = P[X = x_k]@f$. In addition to the
 *         methods specified in the interface
 * @ref umontreal.ssj.probdist.Distribution, a method that returns the probability @f$p_k@f$ is
 *      supplied.
 *
 *      <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_discretereal
 */
public class DiscreteDistribution implements Distribution {
	/*
	 * For better precision in the tails, we keep the cumulative probabilities (F) in cdf[x] for x
	 * <= xmed (i.e. cdf[x] is the sum off all the probabi- lities pr[i] for i <= x), and the
	 * complementary cumulative probabilities (1 - F) in cdf[x] for x > xmed (i.e. cdf[x] is the sum
	 * off all the probabilities pr[i] for i >= x).
	 */

	protected double cdf[] = null;    // cumulative probabilities
	protected double pr[] = null;     // probability terms or mass distribution
	protected int xmin = 0;           // pr[x] = 0 for x < xmin
	protected int xmax = 0;           // pr[x] = 0 for x > xmax
	protected int xmed = 0;           // cdf[x] = F(x) for x <= xmed, and
	                                  // cdf[x] = bar_F(x) for x > xmed
	protected int nVal;               // number of different values
	protected double sortedVal[];
	protected double supportA = Double.NEGATIVE_INFINITY;
	protected double supportB = Double.POSITIVE_INFINITY;

	protected DiscreteDistribution() {
	}
	// Default constructor called by subclasses such as 'EmpiricalDist'

	/**
	 * Constructs a discrete distribution over the @f$n@f$ values contained in array `values`, with
	 * probabilities given in array `prob`. Both arrays must have at least @f$n@f$ elements, the
	 * probabilities must sum to 1, and the values are assumed to be sorted by increasing order.
	 */
	public DiscreteDistribution(double[] values, double[] prob, int n) {
		init(n, values, prob);
	}

	/**
	 * Similar to {@link #DiscreteDistribution(double[],double[],int) DiscreteDistribution(double[],
	 * double[], int)}.
	 */
	public DiscreteDistribution(int[] values, double[] prob, int n) {
		double[] A = new double[n];
		for (int i = 0; i < n; i++)
			A[i] = values[i];
		init(n, A, prob);
	}

	private void init(int n, double[] val, double[] prob) {
		int no = val.length;
		int np = prob.length;
		if (n <= 0)
			throw new IllegalArgumentException("n <= 0");
		if (no < n || np < n)
			throw new IllegalArgumentException("Size of arrays 'values' or 'prob' less than 'n'");

		nVal = n;
		pr = prob;

		// cdf
		sortedVal = new double[nVal];
		System.arraycopy(val, 0, sortedVal, 0, nVal);

		supportA = sortedVal[0];
		supportB = sortedVal[nVal - 1];
		xmin = 0;
		xmax = nVal - 1;

		/*
		 * Compute the cumulative probabilities until F >= 0.5, and keep them in the lower part of
		 * cdf
		 */
		cdf = new double[nVal];
		cdf[0] = pr[0];
		int i = 0;
		while (i < xmax && cdf[i] < 0.5) {
			i++;
			cdf[i] = pr[i] + cdf[i - 1];
		}
		// This is the boundary between F and barF in the CDF
		xmed = i;

		/*
		 * Compute the cumulative probabilities of the complementary distribution and keep them in
		 * the upper part of cdf.
		 */
		cdf[nVal - 1] = pr[nVal - 1];
		i = nVal - 2;
		while (i > xmed) {
			cdf[i] = pr[i] + cdf[i + 1];
			i--;
		}
	}

	/**
	 * @param x
	 *            value at which the cdf is evaluated
	 * @return the cdf evaluated at `x`
	 */
	public double cdf(double x) {
		if (x < sortedVal[0])
			return 0.0;
		if (x >= sortedVal[nVal - 1])
			return 1.0;
		if ((xmax == xmed) || (x < sortedVal[xmed + 1])) {
			for (int i = 0; i <= xmed; i++)
				if (x >= sortedVal[i] && x < sortedVal[i + 1])
					return cdf[i];
		} else {
			for (int i = xmed + 1; i < nVal - 1; i++)
				if (x >= sortedVal[i] && x < sortedVal[i + 1])
					return 1.0 - cdf[i + 1];
		}
		throw new IllegalStateException();
	}

	/**
	 * @param x
	 *            value at which the complementary distribution function is evaluated
	 * @return the complementary distribution function evaluated at `x`
	 */
	public double barF(double x) {
		if (x <= sortedVal[0])
			return 1.0;
		if (x > sortedVal[nVal - 1])
			return 0.0;
		if ((xmax == xmed) || (x <= sortedVal[xmed + 1])) {
			for (int i = 0; i <= xmed; i++)
				if (x > sortedVal[i] && x <= sortedVal[i + 1])
					return 1.0 - cdf[i];
		} else {
			for (int i = xmed + 1; i < nVal - 1; i++)
				if (x > sortedVal[i] && x <= sortedVal[i + 1])
					return cdf[i + 1];
		}
		throw new IllegalStateException();
	}

	/**
	 * @param u
	 *            value in the interval @f$(0,1)@f$ for which the inverse distribution function is
	 *            evaluated
	 * @return the inverse distribution function evaluated at `u`
	 *
	 * @exception IllegalArgumentException
	 *                if @f$u@f$ is not in the interval @f$(0,1)@f$
	 * @exception ArithmeticException
	 *                if the inverse cannot be computed, for example if it would give infinity in a
	 *                theoretical context
	 *
	 */
	public double inverseF(double u) {
		int i, j, k;

		if (u < 0.0 || u > 1.0)
			throw new IllegalArgumentException("u not in [0,1]");
		if (u <= 0.0)
			return supportA;
		if (u >= 1.0)
			return supportB;

		// Remember: the upper part of cdf contains the complementary distribu-
		// tion for xmed < s <= xmax, and the lower part of cdf the
		// distribution for xmin <= s <= xmed

		if (u <= cdf[xmed - xmin]) {
			// In the lower part of cdf
			if (u <= cdf[0])
				return sortedVal[xmin];
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
				return sortedVal[xmax];

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

		return sortedVal[i + xmin];
	}

	/**
	 * Computes the mean @f$E[X] = \sum_i^{} p_i x_i@f$ of the distribution.
	 */
	public double getMean() {
		double mean = 0.0;
		for (int i = 0; i < nVal; i++)
			mean += sortedVal[i] * pr[i];
		return mean;
	}

	/**
	 * Computes the variance @f$\mbox{Var}[X] = \sum_i^{} p_i (x_i - E[X])^2@f$ of the distribution.
	 */
	public double getVariance() {
		double mean = getMean();
		double variance = 0.0;
		for (int i = 0; i < nVal; i++)
			variance += (sortedVal[i] - mean) * (sortedVal[i] - mean) * pr[i];
		return (variance);
	}

	/**
	 * Computes the standard deviation of the distribution.
	 */
	public double getStandardDeviation() {
		return Math.sqrt(getVariance());
	}

	/**
	 * Returns a table containing the parameters of the current distribution. This table is built in
	 * regular order, according to constructor `DiscreteDistribution(double[] params)` order.
	 */
	public double[] getParams() {
		double[] retour = new double[1 + nVal * 2];
		double sum = 0;
		retour[0] = nVal;
		System.arraycopy(sortedVal, 0, retour, 1, nVal);
		for (int i = 0; i < nVal - 1; i++) {
			retour[nVal + 1 + i] = cdf[i] - sum;
			sum = cdf[i];
		}
		retour[2 * nVal] = 1.0 - sum;

		return retour;
	}

	/**
	 * Returns the number of possible values @f$x_i@f$.
	 */
	public int getN() {
		return nVal;
	}

	/**
	 * Returns @f$p_i@f$, the probability of the @f$i@f$-th value, for
	 * 
	 * @f$0\le i<n@f$.
	 * @param i
	 *            value number, @f$0\le i < n@f$
	 * @return the probability of value `i`
	 */
	public double prob(int i) {
		if (i < 0 || i >= nVal)
			return 0.;
		return pr[i];
	}

	/**
	 * Returns the @f$i@f$-th value @f$x_i@f$, for @f$0\le i<n@f$.
	 */
	public double getValue(int i) {
		return sortedVal[i];
	}

	/**
	 * Returns the lower limit @f$x_0@f$ of the support of the distribution.
	 * 
	 * @return @f$x@f$ lower limit of support
	 */
	public double getXinf() {
		return supportA;
	}

	/**
	 * Returns the upper limit @f$x_{n-1}@f$ of the support of the distribution.
	 * 
	 * @return @f$x@f$ upper limit of support
	 */
	public double getXsup() {
		return supportB;
	}

	/**
	 * Returns a `String` containing information about the current distribution.
	 */
	public String toString() {
		StringBuilder sb = new StringBuilder();
		Formatter formatter = new Formatter(sb, Locale.US);
		formatter.format("%s%n", getClass().getSimpleName());
		formatter.format("%s :      %s%n", "value", "cdf");
		for (int i = 0; i < nVal - 1; i++)
			formatter.format("%f : %f%n", sortedVal[i], cdf[i]);
		formatter.format("%f : %f%n", sortedVal[nVal - 1], 1.0);
		return sb.toString();
	}

}