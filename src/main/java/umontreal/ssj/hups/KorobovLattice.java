/*
 * Class:        KorobovLattice
 * Description:  Korobov lattice
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
package umontreal.ssj.hups;

import umontreal.ssj.util.PrintfFormat;

/**
 * This class implements a *Korobov lattice*, which represents the same point set as in class @ref
 * LCGPointSet, but implemented differently. The parameters are the modulus @f$n@f$ and the
 * multiplier @f$a@f$, for an arbitrary integer @f$1 \le a < n@f$. When @f$a@f$ is outside the
 * interval @f$[1,n)@f$, then we replace @f$a@f$ by (@f$a \bmod n@f$) in all calculations. The
 * number of points is @f$n@f$, their dimension is @f$s@f$, and they are defined by @f[ \mathbf{u}_i
 * = (i/n)(1, a, a^2, …, a^{s-1}) \bmod1 @f] for @f$i=0,…,n-1@f$.
 *
 * It is also possible to build a "shifted" Korobov lattice for which the first @f$t@f$ coordinates
 * are skipped. The @f$s@f$-dimensionnal points are then defined as @f[ \mathbf{u}_i = (i/n)(a^t,
 * a^{t+1}, a^{t+2}, …, a^{t+s-1}) \bmod1 @f] for @f$i=0,…,n-1@f$ and fixed @f$t@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class KorobovLattice extends Rank1Lattice {
	protected int genA;            // multiplier a
	private int genT;              // shift t

	private void initN(int n, int t) {
		int a = (genA % n) + (genA < 0 ? n : 0);
		genT = t;
		long[] B = new long[dim];
		B[0] = 1;
		int j;
		for (j = 0; j < t; j++)
			B[0] *= a;
		v[0] = B[0] * normFactor;
		for (j = 1; j < dim; j++) {
			B[j] = (a * B[j - 1]) % n;
			v[j] = normFactor * B[j];
		}
	}

	/**
	 * Instantiates a Korobov lattice point set with modulus @f$n@f$ and multiplier @f$a@f$ in
	 * dimension @f$s@f$.
	 */
	public KorobovLattice(int n, int a, int s) {
		super(n, null, 0);
		genA = a;
		dim = s;
		v = new double[s];
		initN(n, 0);
	}

	/**
	 * Instantiates a shifted Korobov lattice point set with modulus
	 * @f$n@f$ and multiplier @f$a@f$ in dimension @f$s@f$. The first
	 * @f$t@f$ coordinates of a standard Korobov lattice are dropped as described above. The
	 *         case @f$t=0@f$ corresponds to the standard Korobov lattice.
	 */
	public KorobovLattice(int n, int a, int s, int t) {
		super(n, null, 0);
		if (t < 0)
			throw new IllegalArgumentException("KorobovLattice: must have t >= 0");
		dim = s;
		genA = a;
		v = new double[s];
		initN(n, t);
	}

	/**
	 * Resets the number of points of the lattice to @f$n@f$. The values of
	 * @f$s@f$, @f$a@f$ and @f$t@f$ are unchanged.
	 */
	public void setNumPoints(int n) {
		initN(n, genT);
	}

	/**
	 * Returns the multiplier @f$a@f$ of the lattice. (The original one before it is reset to @f$a
	 * \bmod n@f$).
	 */
	public int getA() {
		return genA;
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("KorobovLattice:" + PrintfFormat.NEWLINE);
		sb.append("Multiplier a: " + genA + PrintfFormat.NEWLINE);
		sb.append(super.toString());
		return sb.toString();
	}
}