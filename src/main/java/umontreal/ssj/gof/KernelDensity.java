/*
 * Class:        KernelDensity
 * Description:  Kernel density estimators
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
package umontreal.ssj.gof;

import umontreal.ssj.probdist.*;
// import umontreal.ssj.randvar.KernelDensityGen;

/**
 * This static class provides methods to compute a kernel density estimator from a set of @f$n@f$
 * individual observations @f$x_0, …, x_{n-1}@f$, which define an empirical distribution. The
 * methods return the estimated density either at a given point or on a grid or set of @f$m@f$
 * predefined points. Some details on how a kernel density estimator is defined and how to select
 * the kernel and the bandwidth @f$h@f$ are given in the documentation of class
 * 
 * @ref umontreal.ssj.randvar.KernelDensityGen in package `randvar`.
 *
 *      <div class="SSJ-bigskip"></div>
 */
public class KernelDensity {

	/**
	 * Given an empirical distribution `dist`, this method computes and returns the kernel density
	 * estimate at a single point @f$y@f$. The kernel density is `kern`, and the bandwidth
	 * is @f$h@f$.
	 */
	private static double estimate(EmpiricalDist dist, ContinuousDistribution kern, double h,
	        double y) {
		// Computes and returns the kernel density estimate at $y$, where the
		// kernel is the density kern.density(x), and the bandwidth is $h$.
		double z;
		double a = kern.getXinf();       // lower limit of density
		double b = kern.getXsup();       // upper limit of density
		double sum = 0;
		int n = dist.getN();
		for (int i = 0; i < n; i++) {
			z = (y - dist.getObs(i)) / h;
			if ((z >= a) && (z <= b))
				sum += kern.density(z);
		}
		sum /= (h * n);
		return sum;
	}

	/**
	 * Given the empirical distribution `dist`, this method computes the kernel density estimate at
	 * each of the @f$m@f$ evaluation points <tt>evalPoints[</tt>@f$j@f$<tt>]</tt>, @f$j= 0, 1, …, (m-1)@f$,
	 * where @f$m@f$ is the length of `EvalPoints`, the kernel is `kern.density(x)`, and the bandwidth
	 * is @f$h@f$. Returns the estimates as an array of @f$m@f$ values.
	 * Note that the evaluation points do not have to cover the entire support of the density;
	 * they can cover only a small a small interval.
	 * 
	 * One way to choose @f$h@f$ is via the method
	 * {@link umontreal.ssj.randvar.KernelDensityGen.getBaseBandwidth(EmpiricalDist)
	 * getBaseBandwidth(dist)} in package `randvar`.
	 */
	public static double[] computeDensity(EmpiricalDist dist, ContinuousDistribution kern, double h,
	        double[] Y) {
		int m = Y.length;
		double[] u = new double[m];
		for (int j = 0; j < m; j++)
			u[j] = estimate(dist, kern, h, Y[j]);
		return u;
	}

	/**
	 * Similar to @f$computeDensity (dist, kern, h, evalPoints)@f$, but much more efficient for very
	 * large n. We assume that the kernel is unimodal (increasing, then decreasing). For each
	 * evaluation point j, in the sum that defines the density estimator, we add only the terms that
	 * contribute more than epsilon0 to the sum. We stop adding as soon as the new term is <
	 * epsilon0. One can try epsilon0 = E-10 for example. We have also replaced dist by a sorted
	 * array @f$double[] data@f$ that contains the raw data.
	 * 
	 */
	public void evalDensity(double data[], ContinuousDistribution kern, double h,
	        double[] evalPoints, double[] density, double epsilon0) {
		int m = evalPoints.length;
		// density = new double[m]; // Maybe not needed, but perhaps safer.
		int n = data.length;
		double invhn = 1.0 / (h * n);
		double invh = 1.0 / h;
		double y;
		double sum = 0.0;
		double term;    // A term to be added to the sum that defines the density estimate.
		int imin = 0;   // We know that the terms for i < imin do not contribute significantly.
		for (int j = 0; j < m; j++) {  // Evaluation points are indexed by j.
			y = evalPoints[j];
			term = kern.density((y - data[imin]) * invh);
			while ((term < epsilon0) & (imin < n - 1) && (data[imin] < y))
				term = kern.density((y - data[++imin]) * invh);
			sum = term;   // The first significant term.
			for (int i = imin + 1; (i < n) && ((term > epsilon0) || (data[i] < y)); i++)
				// Data indexed by i.
				sum += (term = kern.density((y - data[i]) * invh));
			density[j] = sum * invhn;
		}
	}

}