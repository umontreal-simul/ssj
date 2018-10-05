/*
 * Class:        DigitalNet
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001--2018  Pierre L'Ecuyer and Universite de Montreal
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
import umontreal.ssj.rng.*;

/**
 * This class provides the basic structures for storing and manipulating <em>linear digital nets in
 * base @f$b@f$</em>, for an arbitrary base @f$b\ge2@f$.
 * 
 * We recall that a net contains @f$n = b^k@f$ points in @f$s@f$ dimensions, where the @f$i@f$th
 * point @f$\mathbf{u}_i@f$, for @f$i=0,…,b^k-1@f$, is defined as follows: @f{align*}{ i & =
 * \sum_{\ell=0}^{k-1} a_{i,\ell} b^{\ell}, \\ \begin{pmatrix} u_{i,j,1} \\ u_{i,j,2} \\ \vdots
 * \end{pmatrix} & = \mathbf{C}_j \begin{pmatrix} a_{i,0} \\ a_{i,1} \\ \vdots \\ a_{i,k-1}
 * \end{pmatrix} , \\ u_{i,j} & = \sum_{\ell=1}^{\infty}u_{i,j,\ell} b^{-\ell}, \\ \mathbf{u}_i & =
 * (u_{i,0},\dots,u_{i,s-1}). @f} In our implementation, the matrices @f$\mathbf{C}_j@f$
 * are @f$r\times k@f$, so the expansion of @f$u_{i,j}@f$ is truncated to its first @f$r@f$ terms.
 * The points are stored implicitly by storing the generator matrices @f$\mathbf{C}_j@f$ in a large
 * two-dimensional array of integers, with @f$srk@f$ elements. For general @f$b@f$, the
 * element @f$(l,c)@f$ of @f$\mathbf{C}_j@f$ (counting elements from 0) is stored at
 * position @f$[jk+c][l]@f$ in this array.
 *
 * To enumerate the points, one should avoid using the method {@link #getCoordinate()
 * getCoordinate(i, j)} for arbitrary values of `i` and `j`, because this is much slower than using
 * a @ref PointSetIterator to access successive coordinates. By default, the iterator enumerates the
 * points @f$\mathbf{u}_i@f$ using a Gray code technique as proposed in @cite rANT79a, @cite
 * rTEZ95a, and also described in @cite fGLA04a, @cite vHON03a). With this technique,
 * the @f$b@f$-ary representation of @f$i@f$, @f$\mathbf{a}_i = (a_{i,0}, …, a_{i,k-1})@f$, is
 * replaced in Equation ({@link REF_hups_overview_eq_digital_Cj digital-Cj}) by a Gray code
 * representation of @f$i@f$, @f$\mathbf{g}_i = (g_{i,0}, …, g_{i,k-1})@f$. The Gray
 * code @f$\mathbf{g}_i@f$ used here is defined by @f$g_{i,k-1} = a_{i,k-1}@f$ and @f$g_{i,\ell} =
 * (a_{i,\ell} - a_{i,\ell+1}) \bmod b@f$ for @f$\ell= 0,…,k-2@f$. It has the property
 * that @f$\mathbf{g}_i = (g_{i,0}, …, g_{i,k-1})@f$ and @f$\mathbf{g}_{i+1} = (g_{i+1,0}, …,
 * g_{i+1,k-1})@f$ differ only in the position of the smallest index @f$\ell@f$ such
 * that @f$a_{i,\ell} < b - 1@f$, and we have @f$g_{i+1,\ell} = (g_{i,\ell}+1) \bmod b@f$ in that
 * position. This Gray code representation permits a more efficient enumeration of the points by the
 * iterators. It changes the order in which the points @f$\mathbf{u}_i@f$ are enumerated, but the
 * first @f$b^m@f$ points remain the same for every integer @f$m@f$. The @f$i@f$th point of the
 * sequence with the Gray enumeration is the @f$i’@f$th point of the original enumeration,
 * where @f$i’@f$ is the integer whose @f$b@f$-ary representation @f$\mathbf{a}_{i’}@f$ is given by
 * the Gray code @f$\mathbf{g}_i@f$. To enumerate all the points successively, we never need to
 * compute the Gray codes explicitly. It suffices to know the position @f$\ell@f$ of the Gray code
 * digit that changes at each step, and this can be found quickly from the @f$b@f$-ary
 * representation @f$\mathbf{a}_i@f$. The digits of each coordinate @f$j@f$ of the current point can
 * be updated by adding column @f$\ell@f$ of the generator matrix @f$\mathbf{C}_j@f$ to the old
 * digits, modulo @f$b@f$.
 *
 * Digital nets can be randomized in various ways @cite mMAT99a, @cite rFAU02a, @cite vLEC02a, @cite
 * vOWE03a. Several types of randomizations specialized for nets are implemented directly in this
 * class. A simple but important randomization is the *random digital shift* in base @f$b@f$,
 * defined as follows: replace each digit @f$u_{i,j,\ell}@f$ in (
 * {@link REF_hups_overview_eq_digital_uij digital-uij} ) by @f$(u_{i,j,\ell} + d_{j,\ell}) \bmod
 * b@f$, where the @f$d_{j,\ell}@f$’s are i.i.d. uniform over @f$\{0,\dots,b-1\}@f$. This is
 * equivalent to applying a single random shift to all the points in a formal series representation
 * of their coordinates @cite vLEC02a, @cite vLEM03a. In practice, the digital shift is truncated
 * to @f$w@f$ digits, for some integer @f$w\ge r@f$. Applying a digital shift does not change the
 * equidistribution and @f$(t,m,s)@f$-net properties of a point set @cite vHON03a, @cite
 * vLEC99a, @cite vLEM03a. Moreover, with the random shift, each point has the uniform distribution
 * over the unit hypercube (but the points are not independent, of course).
 *
 * A second class of randomizations specialized for digital nets are the *linear matrix
 * scrambles* @cite mMAT99a, @cite rFAU02a, @cite vHON03a, @cite vOWE03a, which multiply the
 * matrices @f$\mathbf{C}_j@f$ by a random invertible matrix @f$\mathbf{M}_j@f$, modulo @f$b@f$.
 * There are several variants, depending on how @f$\mathbf{M}_j@f$ is generated, and on
 * whether @f$\mathbf{C}_j@f$ is multiplied on the left or on the right. In our implementation, the
 * linear matrix scrambles are incorporated directly into the matrices @f$\mathbf{C}_j@f$ (as
 * in @cite vHON03a), so they do not slow down the enumeration of points. Methods are available for
 * applying linear matrix scrambles and for removing these randomizations. These methods generate
 * the appropriate random numbers and make the corresponding changes to the @f$\mathbf{C}_j@f$’s. A
 * copy of the original \f$\mathbf{C}_j\f$’s is maintained, so the point set can be returned to its
 * original unscrambled state at any time. When a new linear matrix scramble is applied, it is
 * always applied to the *original* generator matrices. The method #resetGeneratorMatrices removes
 * the current matrix scramble by resetting the generator matrices to their original state. On the
 * other hand, the method #eraseOriginalGeneratorMatrices replaces the original generator matrices
 * by the current ones, making the changes permanent. This could be useful if one wishes to apply
 * two or more linear matrix scrambles on top of each other and not retain the original matrices.
 *
 * With the linear matrix scrambles alone, the randomized points do not have the uniform
 * distribution over the unit cube. For this reason, they are usually combined with a random digital
 * shift; this combination is called an *affine matrix scramble* @cite vOWE03a. These two
 * randomizations are applied via separate methods. The linear matrix scrambles are incorporated
 * into the matrices @f$\mathbf{C}_j@f$ whereas the digital random shift is stored and applied
 * separately, independently of the other scrambles.
 *
 * Applying a digital shift or a linear matrix scramble to a digital net invalidates all current
 * iterators for the current point, because each iterator uses a *cached* copy of the current point,
 * which is updated only when the current point index of that iterator changes, and the update also
 * depends on the cached copy of the previous point. After applying any kind of scrambling or
 * randomization that affects the `DigitalNet` object, the iterators must be reinitialized to the
 * *initial point* by invoking `PointSetIterator.resetCurPointIndex` or re-instantiated by the
 * `iterator` method (this is not done automatically).
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class DigitalNet extends PointSet {

	// Variables to be initialized by the constructor subclasses.
	protected int b = 0;         // Base.
	protected int numCols = 0;   // The number of columns in each C_j. (= k)
	protected int numRows = 0;   // The number of rows in each C_j. (= r)
	protected int outDigits = 0; // Number of output digits (= w)
	private int[][] originalMat; // Original gen. matrices without randomizat.
	protected int[][] genMat;    // The current generator matrices.
	                             // genMat[j*numCols+c][l] contains column c
	                             // and row l of C_j.
	protected int[][] digitalShift; // The digital shift, initially zero (null).
	                                // Entry [j][l] is for dimension j, digit l,
	                                // for 0 <= l < outDigits.
	protected double normFactor; // To convert output to (0,1); 1/b^outDigits.
	protected double[] factor;   // Lookup table in ascending order: factor[i]
	                             // = 1/b^{i+1} for 0 <= i < outDigits.

	// primes gives the first index in array FaureFactor
	// for the prime p. If primes[i] = p, then
	// FaureFactor[p][j] contains the Faure ordered factors of base p.
	private int[] primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
	        67 };

	// Factors on the diagonal corresponding to base b = prime[i] ordered by
	// increasing Bounds.
	private int[][] FaureFactor = { { 1 }, { 1, 2 }, { 2, 3, 1, 4 }, { 2, 3, 4, 5, 1, 6 },
	        { 3, 4, 7, 8, 2, 5, 6, 9, 1, 10 }, { 5, 8, 3, 4, 9, 10, 2, 6, 7, 11, 1, 12 },
	        { 5, 7, 10, 12, 3, 6, 11, 14, 4, 13, 2, 8, 9, 15, 1, 16 },
	        { 7, 8, 11, 12, 4, 5, 14, 15, 3, 6, 13, 16, 2, 9, 10, 17, 1, 18 },
	        { 5, 9, 14, 18, 7, 10, 13, 16, 4, 6, 17, 19, 3, 8, 15, 20, 2, 11, 12, 21, 1, 22 },
	        { 8, 11, 18, 21, 12, 17, 9, 13, 16, 20, 5, 6, 23, 24, 4, 7, 22, 25, 3, 10, 19, 26, 2,
	                14, 15, 27, 1, 28 },
	        { 12, 13, 18, 19, 11, 14, 17, 20, 7, 9, 22, 24, 4, 8, 23, 27, 5, 6, 25, 26, 3, 10, 21,
	                28, 2, 15, 16, 29, 1, 30 },
	        { 8, 14, 23, 29, 10, 11, 26, 27, 13, 17, 20, 24, 7, 16, 21, 30, 5, 15, 22, 32, 6, 31, 4,
	                9, 28, 33, 3, 12, 25, 34, 2, 18, 19, 35, 1, 36 },
	        { 16, 18, 23, 25, 11, 15, 26, 30, 12, 17, 24, 29, 9, 32, 13, 19, 22, 28, 6, 7, 34, 35,
	                5, 8, 33, 36, 4, 10, 31, 37, 3, 14, 27, 38, 2, 20, 21, 39, 1, 40 },
	        { 12, 18, 25, 31, 9, 19, 24, 34, 8, 16, 27, 35, 10, 13, 30, 33, 15, 20, 23, 28, 5, 17,
	                26, 38, 6, 7, 36, 37, 4, 11, 32, 39, 3, 14, 29, 40, 2, 21, 22, 41, 1, 42 },
	        { 13, 18, 29, 34, 11, 17, 30, 36, 10, 14, 33, 37, 7, 20, 27, 40, 9, 21, 26, 38, 15, 22,
	                25, 32, 6, 8, 39, 41, 5, 19, 28, 42, 4, 12, 35, 43, 3, 16, 31, 44, 2, 23, 24,
	                45, 1, 46 },
	        { 14, 19, 34, 39, 23, 30, 12, 22, 31, 41, 8, 11, 20, 24, 29, 33, 42, 45, 10, 16, 37, 43,
	                7, 15, 38, 46, 17, 25, 28, 36, 5, 21, 32, 48, 6, 9, 44, 47, 4, 13, 40, 49, 3,
	                18, 35, 50, 2, 26, 27, 51, 1, 52 },
	        { 25, 26, 33, 34, 18, 23, 36, 41, 14, 21, 38, 45, 24, 27, 32, 35, 11, 16, 43, 48, 9, 13,
	                46, 50, 8, 22, 37, 51, 7, 17, 42, 52, 19, 28, 31, 40, 6, 10, 49, 53, 5, 12, 47,
	                54, 4, 15, 44, 55, 3, 20, 39, 56, 2, 29, 30, 57, 1, 58 },
	        { 22, 25, 36, 39, 17, 18, 43, 44, 24, 28, 33, 37, 13, 14, 47, 48, 16, 19, 42, 45, 9, 27,
	                34, 52, 8, 23, 38, 53, 11, 50, 7, 26, 35, 54, 21, 29, 32, 40, 6, 10, 51, 55, 5,
	                12, 49, 56, 4, 15, 46, 57, 3, 20, 41, 58, 2, 30, 31, 59, 1, 60 },
	        { 18, 26, 41, 49, 14, 24, 43, 53, 12, 28, 39, 55, 29, 30, 37, 38, 10, 20, 47, 57, 16,
	                21, 46, 51, 8, 25, 42, 59, 13, 31, 36, 54, 9, 15, 52, 58, 7, 19, 48, 60, 23, 32,
	                35, 44, 5, 27, 40, 62, 6, 11, 56, 61, 4, 17, 50, 63, 3, 22, 45, 64, 2, 33, 34,
	                65, 1, 66 } };

	/**
	 * Empty constructor.
	 */
	public DigitalNet() {
	}

	/**
	 * Returns @f$u_{i',j}@f$, the coordinate @f$j@f$ of point @f$i'@f$, where @f$i'@f$ is the Gray
	 * code for @f$i@f$.
	 * 
	 * @param i
	 *            point index, to be transformed to a Gray code
	 * @param j
	 *            coordinate index
	 * @return the value of @f$u_{i,j}@f$
	 */
	public double getCoordinate(int i, int j) {
		// convert i to Gray representation, put digits in gdigit[].
		int l, c, sum;
		int[] bdigit = new int[numCols];
		int[] gdigit = new int[numCols];
		int idigits = intToDigitsGray(b, i, numCols, bdigit, gdigit);
		double result = 0;
		if (digitalShift != null && dimShift < j)
			addRandomShift(dimShift, j, shiftStream);  // Extends the random shift up to j if
			                                           // needed.
		for (l = 0; l < outDigits; l++) {
			if (digitalShift == null)
				sum = 0;
			else
				sum = digitalShift[j][l];
			if (l < numRows)
				for (c = 0; c < idigits; c++)
					sum += genMat[j * numCols + c][l] * gdigit[c];
			result += (sum % b) * factor[l];
		}
		if (digitalShift != null)
			result += EpsilonHalf;
		return result;
	}

	/**
	 * @return an iterator to this point set, using the Gray code enumeration.
	 */
	public PointSetIterator iterator() {
		return new DigitalNetIterator();
	}

	/**
	 * Returns @f$u_{i,j}@f$, the coordinate @f$j@f$ of point @f$i@f$, the points being enumerated
	 * in the standard order (no Gray code).
	 * 
	 * @param i
	 *            point index
	 * @param j
	 *            coordinate index
	 * @return the value of @f$u_{i,j}@f$
	 */
	public double getCoordinateNoGray(int i, int j) {
		// convert i to b-ary representation, put digits in bdigit[].
		int l, c, sum;
		int[] bdigit = new int[numCols];
		int idigits = 0;
		for (c = 0; i > 0; c++) {
			idigits++;
			bdigit[c] = i % b;
			i = i / b;
		}
		if (digitalShift != null && dimShift < j)
			addRandomShift(dimShift, j, shiftStream); // Extends the random shift up to j if needed.
		double result = 0;
		for (l = 0; l < outDigits; l++) {
			if (digitalShift == null)
				sum = 0;
			else
				sum = digitalShift[j][l];
			if (l < numRows)
				for (c = 0; c < idigits; c++)
					sum += genMat[j * numCols + c][l] * bdigit[c];
			result += (sum % b) * factor[l];
		}
		if (digitalShift != null)
			result += EpsilonHalf;
		return result;
	}

	/**
	 * Returns an iterator that does not use the Gray code. With this one, the points will be
	 * enumerated in the order of their first coordinate before randomization.
	 */
	public PointSetIterator iteratorNoGray() {
		return new DigitalNetIteratorNoGray();
	}

	/**
	 * Generates a random digital shift for coordinates @f$j@f$ from `d1` to `d2-1`, using `stream`
	 * to generate the random numbers. The dimension of the current shift is reset to `d2` and the
	 * current `streamShift` is set to `stream`. This shift vector @f$(d_{j,0},…,d_{j,k-1})@f$ is
	 * generated uniformly over @f$\{0,\dots,b-1\}^k@f$ for each coordinate. This shift vector will
	 * be added modulo @f$b@f$ to the digits of all the points by any iterator on this point set.
	 * After adding a digital shift, all iterators must be reconstructed or reset to zero.
	 * 
	 * @param stream
	 *            random number stream used to generate the uniforms
	 */
	public void addRandomShift(int d1, int d2, RandomStream stream) {
		if (null == stream)
			throw new IllegalArgumentException(
			        PrintfFormat.NEWLINE + "   Calling addRandomShift with null stream");
		if (0 == d2)
			d2 = Math.max(1, dim);
		if (digitalShift == null) {
			digitalShift = new int[d2][outDigits];
			capacityShift = d2;
		} else if (d2 > capacityShift) {
			int d3 = Math.max(4, capacityShift);
			while (d2 > d3)
				d3 *= 2;
			int[][] temp = new int[d3][outDigits];
			capacityShift = d3;
			for (int i = 0; i < d1; i++)
				for (int j = 0; j < outDigits; j++)
					temp[i][j] = digitalShift[i][j];
			digitalShift = temp;
		}
		for (int i = d1; i < d2; i++)
			for (int j = 0; j < outDigits; j++)
				digitalShift[i][j] = stream.nextInt(0, b - 1);
		dimShift = d2;
		shiftStream = stream;
	}

	/**
	 * Same as {@link #addRandomShift() addRandomShift(0, dim, stream)}, where `dim` is the
	 * dimension of the digital net.
	 * 
	 * @param stream
	 *            random number stream used to generate the uniforms
	 */
	public void addRandomShift(RandomStream stream) {
		addRandomShift(0, dim, stream);
	}

	/**
	 * Erases the current digital random shift, if any.
	 */
	public void clearRandomShift() {
		super.clearRandomShift();
		digitalShift = null;
	}

	/**
	 * Formats a string that contains information on this digital net.
	 * 
	 * @return string representation of basic information on this digital net
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer(100);
		if (b > 0) {
			sb.append(", base = ");
			sb.append(b);
		}
		sb.append(", Num cols = ");
		sb.append(numCols);
		sb.append(", Num rows = ");
		sb.append(numRows);
		sb.append(", outDigits = ");
		sb.append(outDigits);
		return sb.toString();
	}

	// Print matrices M for dimensions 0 to N-1.
	protected void printMat(int N, int[][][] A, String name) {
		for (int i = 0; i < N; i++) {
			System.out.println("-------------------------------------" + PrintfFormat.NEWLINE + name
			        + "   dim = " + i);
			int l, c;   // row l, column c, dimension i for A[i][l][c].
			for (l = 0; l < numRows; l++) {
				for (c = 0; c < numCols; c++) {
					System.out.print(A[i][l][c] + "  ");
				}
				System.out.println("");
			}
		}
		System.out.println("");
	}

	// Print matrix M
	protected void printMat0(int[][] A, String name) {
		System.out.println("-------------------------------------" + PrintfFormat.NEWLINE + name);
		int l, c;   // row l, column c for A[l][c].
		for (l = 0; l < numCols; l++) {
			for (c = 0; c < numCols; c++) {
				System.out.print(A[l][c] + "  ");
			}
			System.out.println("");
		}
		System.out.println("");
	}

	// Left-multiplies lower-triangular matrix Mj by original C_j,
	// where original C_j is in originalMat and result is in genMat.
	// This implementation is safe only if (numRows*(b-1)^2) is a valid int.
	private void leftMultiplyMat(int j, int[][] Mj) {
		int l, c, i, sum;   // Dimension j, row l, column c for new C_j.
		for (l = 0; l < numRows; l++) {
			for (c = 0; c < numCols; c++) {
				// Multiply row l of M_j by column c of C_j.
				sum = 0;
				for (i = 0; i <= l; i++)
					sum += Mj[l][i] * originalMat[j * numCols + c][i];
				genMat[j * numCols + c][l] = sum % b;
			}
		}
	}

	// Left-multiplies diagonal matrix Mj by original C_j,
	// where original C_j is in originalMat and result is in genMat.
	// This implementation is safe only if (numRows*(b-1)^2) is a valid int.
	private void leftMultiplyMatDiag(int j, int[][] Mj) {
		int l, c, sum;   // Dimension j, row l, column c for new C_j.
		for (l = 0; l < numRows; l++) {
			for (c = 0; c < numCols; c++) {
				// Multiply row l of M_j by column c of C_j.
				sum = Mj[l][l] * originalMat[j * numCols + c][l];
				genMat[j * numCols + c][l] = sum % b;
			}
		}
	}

	// Right-multiplies original C_j by upper-triangular matrix Mj,
	// where original C_j is in originalMat and result is in genMat.
	// This implementation is safe only if (numCols*(b-1)^2) is a valid int.
	private void rightMultiplyMat(int j, int[][] Mj) {
		int l, c, i, sum;   // Dimension j, row l, column c for new C_j.
		for (l = 0; l < numRows; l++) {
			for (c = 0; c < numCols; c++) {
				// Multiply row l of C_j by column c of M_j.
				sum = 0;
				for (i = 0; i <= c; i++)
					sum += originalMat[j * numCols + i][l] * Mj[i][c];
				genMat[j * numCols + c][l] = sum % b;
			}
		}
	}

	private int getFaureIndex(String method, int sb, int flag) {
		// Check for errors in ...FaurePermut. Returns the index ib of the
		// base b in primes, i.e. b = primes[ib].
		if (sb >= b)
			throw new IllegalArgumentException(PrintfFormat.NEWLINE + "   sb >= base in " + method);
		if (sb < 1)
			throw new IllegalArgumentException(PrintfFormat.NEWLINE + "   sb = 0 in " + method);
		if ((flag > 2) || (flag < 0))
			throw new IllegalArgumentException(
			        PrintfFormat.NEWLINE + "   lowerFlag not in {0, 1, 2} in " + method);

		// Find index ib of base in array primes
		int ib = 0;
		while ((ib < primes.length) && (primes[ib] < b))
			ib++;
		if (ib >= primes.length)
			throw new IllegalArgumentException("base too large in " + method);
		if (b != primes[ib])
			throw new IllegalArgumentException(
			        "Faure factors are not implemented for this base in " + method);
		return ib;
	}

	/**
	 * Applies a linear scramble by multiplying each @f$\mathbf{C}_j@f$ on the left by a @f$w\times
	 * w@f$ nonsingular lower-triangular matrix
	 * 
	 * @f$\mathbf{M}_j@f$, as suggested by Matoušek @cite mMAT99a&thinsp; and implemented by Hong
	 *                     and Hickernell @cite vHON03a&thinsp;. The diagonal entries of each
	 *                     matrix @f$\mathbf{M}_j@f$ are generated uniformly over @f$\{1,…,b-1\}@f$,
	 *                     the entries below the diagonal are generated uniformly
	 *                     over @f$\{0,…,b-1\}@f$, and all these entries are generated
	 *                     independently. This means that in base @f$b=2@f$, all diagonal elements
	 *                     are equal to 1.
	 * @remark **Richard:** Les matrices de `leftMatrixScramble` sont carrées et triangulaires
	 *         inférieures. PL pense qu’il faut considérer la possibilité de rajouter des lignes à
	 *         ces matrices pour pouvoir randomiser plus les derniers chiffres ou les derniers bits.
	 * @param stream
	 *            random number stream used to generate the randomness
	 */
	public void leftMatrixScramble(RandomStream stream) {
		int j, l, c;  // dimension j, row l, column c.

		// If genMat contains the original gen. matrices, copy to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Constructs the lower-triangular scrambling matrices M_j, r by r.
		int[][][] scrambleMat = new int[dim][numRows][numRows];
		for (j = 0; j < dim; j++) {
			for (l = 0; l < numRows; l++) {
				for (c = 0; c < numRows; c++) {
					if (c == l)                   // No zero on the diagonal.
						scrambleMat[j][l][c] = stream.nextInt(1, b - 1);
					else if (c < l)
						scrambleMat[j][l][c] = stream.nextInt(0, b - 1);
					else
						scrambleMat[j][l][c] = 0;  // Zeros above the diagonal;
				}
			}
		}

		// Multiply M_j by the generator matrix C_j for each j.
		for (j = 0; j < dim; j++)
			leftMultiplyMat(j, scrambleMat[j]);
	}

	/**
	 * Similar to #leftMatrixScramble except that all the off-diagonal elements of
	 * the @f$\mathbf{M}_j@f$ are 0.
	 * 
	 * @param stream
	 *            random number stream used to generate the randomness
	 */
	public void leftMatrixScrambleDiag(RandomStream stream) {
		int j, l, c;  // dimension j, row l, column c.

		// If genMat contains the original gen. matrices, copy to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Constructs the diagonal scrambling matrices M_j, r by r.
		int[][][] scrambleMat = new int[dim][numRows][numRows];
		for (j = 0; j < dim; j++) {
			for (l = 0; l < numRows; l++) {
				for (c = 0; c < numRows; c++) {
					if (c == l)                   // No zero on the diagonal.
						scrambleMat[j][l][c] = stream.nextInt(1, b - 1);
					else
						scrambleMat[j][l][c] = 0;  // Diagonal matrix;
				}
			}
		}

		// Multiply M_j by the generator matrix C_j for each j.
		for (j = 0; j < dim; j++)
			leftMultiplyMatDiag(j, scrambleMat[j]);
	}

	private void LMSFaurePermut(String method, RandomStream stream, int sb, int lowerFlag) {
		/*
		 * If \texttt{lowerFlag = 2}, the off-diagonal elements below the diagonal of $\mathbf{M}_j$
		 * are chosen as in \method{leftMatrixScramble}{}. If \texttt{lowerFlag = 1}, the
		 * off-diagonal elements below the diagonal of $\mathbf{M}_j$ are also chosen from the
		 * restricted set. If \texttt{lowerFlag = 0}, the off-diagonal elements of $\mathbf{M}_j$
		 * are all 0.
		 */
		int ib = getFaureIndex(method, sb, lowerFlag);

		// If genMat contains the original gen. matrices, copy to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Constructs the lower-triangular scrambling matrices M_j, r by r.
		int jb;
		int j, l, c;  // dimension j, row l, column c.
		int[][][] scrambleMat = new int[dim][numRows][numRows];
		for (j = 0; j < dim; j++) {
			for (l = 0; l < numRows; l++) {
				for (c = 0; c < numRows; c++) {
					if (c == l) {
						jb = stream.nextInt(0, sb - 1);
						scrambleMat[j][l][c] = FaureFactor[ib][jb];
					} else if (c < l) {
						if (lowerFlag == 2) {
							scrambleMat[j][l][c] = stream.nextInt(0, b - 1);
						} else if (lowerFlag == 1) {
							jb = stream.nextInt(0, sb - 1);
							scrambleMat[j][l][c] = FaureFactor[ib][jb];
						} else {   // lowerFlag == 0
							scrambleMat[j][l][c] = 0;
						}
					} else
						scrambleMat[j][l][c] = 0;  // Zeros above the diagonal;
				}
			}
		}
		// printMat (dim, scrambleMat, method);

		// Multiply M_j by the generator matrix C_j for each j.
		if (lowerFlag == 0)
			for (j = 0; j < dim; j++)
				leftMultiplyMatDiag(j, scrambleMat[j]);
		else
			for (j = 0; j < dim; j++)
				leftMultiplyMat(j, scrambleMat[j]);
	}

	/**
	 * Similar to #leftMatrixScramble except that the diagonal elements of each
	 * matrix @f$\mathbf{M}_j@f$ are chosen from a restricted set of the best integers as calculated
	 * by Faure @cite rFAU02a&thinsp;. They are generated uniformly over the first `sb` elements of
	 * array @f$F@f$, where @f$F@f$ is made up of a permutation of the integers @f$[1..(b-1)]@f$.
	 * These integers are sorted by increasing order of the upper bounds of the extreme discrepancy
	 * for the given integer.
	 * 
	 * @param stream
	 *            random number stream used to generate the randomness
	 * @param sb
	 *            Only the first @f$sb@f$ elements of @f$F@f$ are used
	 */
	public void leftMatrixScrambleFaurePermut(RandomStream stream, int sb) {
		LMSFaurePermut("leftMatrixScrambleFaurePermut", stream, sb, 2);
	}

	/**
	 * Similar to #leftMatrixScrambleFaurePermut except that all off-diagonal elements are 0.
	 * 
	 * @param stream
	 *            random number stream used to generate the randomness
	 * @param sb
	 *            Only the first @f$sb@f$ elements of @f$F@f$ are used
	 */
	public void leftMatrixScrambleFaurePermutDiag(RandomStream stream, int sb) {
		LMSFaurePermut("leftMatrixScrambleFaurePermutDiag", stream, sb, 0);
	}

	/**
	 * Similar to #leftMatrixScrambleFaurePermut except that the elements under the diagonal are
	 * also chosen from the same restricted set as the diagonal elements.
	 * 
	 * @param stream
	 *            random number stream used to generate the randomness
	 * @param sb
	 *            Only the first @f$sb@f$ elements of @f$F@f$ are used
	 */
	public void leftMatrixScrambleFaurePermutAll(RandomStream stream, int sb) {
		LMSFaurePermut("leftMatrixScrambleFaurePermutAll", stream, sb, 1);
	}

	/**
	 * Applies the @f$i@f$-binomial matrix scramble proposed by Tezuka
	 * 
	 * @cite rTEZ02a&thinsp; (see also @cite vOWE03a&thinsp;). This multiplies
	 *       each @f$\mathbf{C}_j@f$ on the left by a @f$w\times w@f$ nonsingular lower-triangular
	 *       matrix @f$\mathbf{M}_j@f$ as in #leftMatrixScramble, but with the additional constraint
	 *       that all entries on any given diagonal or subdiagonal of @f$\mathbf{M}_j@f$ are
	 *       identical.
	 * @param stream
	 *            random number stream used as generator of the randomness
	 */
	public void iBinomialMatrixScramble(RandomStream stream) {
		int j, l, c;  // dimension j, row l, column c.
		int diag;     // random entry on the diagonal;
		int col1;     // random entries below the diagonal;

		// If genMat is original generator matrices, copy it to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Constructs the lower-triangular scrambling matrices M_j, r by r.
		int[][][] scrambleMat = new int[dim][numRows][numRows];
		for (j = 0; j < dim; j++) {
			diag = stream.nextInt(1, b - 1);
			for (l = 0; l < numRows; l++) {
				// Single random nonzero element on the diagonal.
				scrambleMat[j][l][l] = diag;
				for (c = l + 1; c < numRows; c++)
					scrambleMat[j][l][c] = 0;
			}
			for (l = 1; l < numRows; l++) {
				col1 = stream.nextInt(0, b - 1);
				for (c = 0; l + c < numRows; c++)
					scrambleMat[j][l + c][c] = col1;
			}
		}
		// printMat (dim, scrambleMat, "iBinomialMatrixScramble");
		for (j = 0; j < dim; j++)
			leftMultiplyMat(j, scrambleMat[j]);
	}

	private void iBMSFaurePermut(String method, RandomStream stream, int sb, int lowerFlag) {
		int ib = getFaureIndex(method, sb, lowerFlag);

		// If genMat is original generator matrices, copy it to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Constructs the lower-triangular scrambling matrices M_j, r by r.
		int j, l, c;  // dimension j, row l, column c.
		int diag;     // random entry on the diagonal;
		int col1;     // random entries below the diagonal;
		int jb;
		int[][][] scrambleMat = new int[dim][numRows][numRows];
		for (j = 0; j < dim; j++) {
			jb = stream.nextInt(0, sb - 1);
			diag = FaureFactor[ib][jb];
			for (l = 0; l < numRows; l++) {
				// Single random nonzero element on the diagonal.
				scrambleMat[j][l][l] = diag;
				for (c = l + 1; c < numRows; c++)
					scrambleMat[j][l][c] = 0;
			}
			for (l = 1; l < numRows; l++) {
				if (lowerFlag == 2) {
					col1 = stream.nextInt(0, b - 1);
				} else if (lowerFlag == 1) {
					jb = stream.nextInt(0, sb - 1);
					col1 = FaureFactor[ib][jb];
				} else {   // lowerFlag == 0
					col1 = 0;
				}
				for (c = 0; l + c < numRows; c++)
					scrambleMat[j][l + c][c] = col1;
			}
		}
		// printMat (dim, scrambleMat, method);

		if (lowerFlag > 0)
			for (j = 0; j < dim; j++)
				leftMultiplyMat(j, scrambleMat[j]);
		else
			for (j = 0; j < dim; j++)
				leftMultiplyMatDiag(j, scrambleMat[j]);
	}

	/**
	 * Similar to #iBinomialMatrixScramble except that the diagonal elements of each
	 * matrix @f$\mathbf{M}_j@f$ are chosen as in #leftMatrixScrambleFaurePermut.
	 * 
	 * @param stream
	 *            random number stream used to generate the randomness
	 * @param sb
	 *            Only the first @f$sb@f$ elements of @f$F@f$ are used
	 */
	public void iBinomialMatrixScrambleFaurePermut(RandomStream stream, int sb) {
		iBMSFaurePermut("iBinomialMatrixScrambleFaurePermut", stream, sb, 2);
	}

	/**
	 * Similar to #iBinomialMatrixScrambleFaurePermut except that all the off-diagonal elements are
	 * 0.
	 * 
	 * @param stream
	 *            random number stream used to generate the randomness
	 * @param sb
	 *            Only the first @f$sb@f$ elements of @f$F@f$ are used
	 */
	public void iBinomialMatrixScrambleFaurePermutDiag(RandomStream stream, int sb) {
		iBMSFaurePermut("iBinomialMatrixScrambleFaurePermutDiag", stream, sb, 0);
	}

	/**
	 * Similar to #iBinomialMatrixScrambleFaurePermut except that the elements under the diagonal
	 * are also chosen from the same restricted set as the diagonal elements.
	 * 
	 * @param stream
	 *            random number stream used to generate the randomness
	 * @param sb
	 *            Only the first @f$sb@f$ elements of @f$F@f$ are used
	 */
	public void iBinomialMatrixScrambleFaurePermutAll(RandomStream stream, int sb) {
		iBMSFaurePermut("iBinomialMatrixScrambleFaurePermutAll", stream, sb, 1);
	}

	/**
	 * Applies the *striped matrix scramble* proposed by Owen
	 * 
	 * @cite vOWE03a&thinsp;. It multiplies each @f$\mathbf{C}_j@f$ on the left by a @f$w\times w@f$
	 *       nonsingular lower-triangular matrix
	 * @f$\mathbf{M}_j@f$ as in #leftMatrixScramble, but with the additional constraint that in any
	 *                    column, all entries below the diagonal are equal to the diagonal entry,
	 *                    which is generated randomly over @f$\{1,…,b-1\}@f$. Note that
	 *                    for @f$b=2@f$, the matrices @f$\mathbf{M}_j@f$ become deterministic, with
	 *                    all entries on and below the diagonal equal to 1.
	 * @param stream
	 *            random number stream used as generator of the randomness
	 */
	public void stripedMatrixScramble(RandomStream stream) {
		int j, l, c;  // dimension j, row l, column c.
		int diag;     // random entry on the diagonal;
		// int col1; // random entries below the diagonal;

		// If genMat contains original gener. matrices, copy it to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Constructs the lower-triangular scrambling matrices M_j, r by r.
		int[][][] scrambleMat = new int[dim][numRows][numRows];
		for (j = 0; j < dim; j++) {
			for (c = 0; c < numRows; c++) {
				diag = stream.nextInt(1, b - 1);   // Random entry in this column.
				for (l = 0; l < c; l++)
					scrambleMat[j][l][c] = 0;
				for (l = c; l < numRows; l++)
					scrambleMat[j][l][c] = diag;
			}
		}
		// printMat (dim, scrambleMat, "stripedMatrixScramble");
		for (j = 0; j < dim; j++)
			leftMultiplyMat(j, scrambleMat[j]);
	}

	/**
	 * Similar to #stripedMatrixScramble except that the elements on and under the diagonal of each
	 * matrix @f$\mathbf{M}_j@f$ are chosen as in #leftMatrixScrambleFaurePermut.
	 * 
	 * @param stream
	 *            random number stream used as generator of the randomness
	 * @param sb
	 *            Only the first @f$sb@f$ elements of @f$F@f$ are used
	 */
	public void stripedMatrixScrambleFaurePermutAll(RandomStream stream, int sb) {
		int ib = getFaureIndex("stripedMatrixScrambleFaurePermutAll", sb, 1);

		// If genMat contains original gener. matrices, copy it to originalMat.
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Constructs the lower-triangular scrambling matrices M_j, r by r.
		int j, l, c;  // dimension j, row l, column c.
		int diag;     // random entry on the diagonal;
		// int col1; // random entries below the diagonal;
		int jb;
		int[][][] scrambleMat = new int[dim][numRows][numRows];
		for (j = 0; j < dim; j++) {
			for (c = 0; c < numRows; c++) {
				jb = stream.nextInt(0, sb - 1);
				diag = FaureFactor[ib][jb];  // Random entry in this column.
				for (l = 0; l < c; l++)
					scrambleMat[j][l][c] = 0;
				for (l = c; l < numRows; l++)
					scrambleMat[j][l][c] = diag;
			}
		}
		// printMat (dim, scrambleMat, "stripedMatrixScrambleFaurePermutAll");
		for (j = 0; j < dim; j++)
			leftMultiplyMat(j, scrambleMat[j]);
	}

	/**
	 * Applies a linear scramble by multiplying each @f$\mathbf{C}_j@f$ on the right by a
	 * single @f$k\times k@f$ nonsingular upper-triangular matrix @f$\mathbf{M}@f$, as suggested by
	 * Faure and Tezuka
	 * 
	 * @cite rFAU02a&thinsp; (see also @cite vHON03a&thinsp;). The diagonal entries of the
	 *       matrix @f$\mathbf{M}@f$ are generated uniformly over @f$\{1,…,b-1\}@f$, the entries
	 *       above the diagonal are generated uniformly over @f$\{0,…,b-1\}@f$, and all the entries
	 *       are generated independently. The effect of this scramble is only to change the order in
	 *       which the points are generated. If one computes the average value of a function over
	 *       *all* the points of a given digital net, or over a number of points that is a power of
	 *       the basis, then this scramble makes no difference.
	 * @param stream
	 *            random number stream used as generator of the randomness
	 */
	public void rightMatrixScramble(RandomStream stream) {
		int j, c, l;  // dimension j, row l, column c, of genMat.

		// SaveOriginalMat();
		if (originalMat == null) {
			originalMat = genMat;
			genMat = new int[dim * numCols][numRows];
		}

		// Generate an upper triangular matrix for Faure-Tezuka right-scramble.
		// Entry [l][c] is in row l, col c.
		int[][] scrambleMat = new int[numCols][numCols];
		for (c = 0; c < numCols; c++) {
			for (l = 0; l < c; l++)
				scrambleMat[l][c] = stream.nextInt(0, b - 1);
			scrambleMat[c][c] = stream.nextInt(1, b - 1);
			for (l = c + 1; l < numCols; l++)
				scrambleMat[l][c] = 0;
		}

		// printMat0 (scrambleMat, "rightMatrixScramble");
		// Right-multiply the generator matrices by the scrambling matrix.
		for (j = 0; j < dim; j++)
			rightMultiplyMat(j, scrambleMat);
	}

	/**
	 * Restores the original generator matrices and removes the random shift.
	 */
	public void unrandomize() {
		resetGeneratorMatrices();
		digitalShift = null;
	}

	/**
	 * Restores the original generator matrices. This removes the current linear matrix scrambles.
	 */
	public void resetGeneratorMatrices() {
		if (originalMat != null) {
			genMat = originalMat;
			originalMat = null;
		}
	}

	/**
	 * Erases the original generator matrices and replaces them by the current ones. The current
	 * linear matrix scrambles thus become *permanent*. This is useful if we want to apply several
	 * scrambles in succession to a given digital net.
	 */
	public void eraseOriginalGeneratorMatrices() {
		originalMat = null;
	}

	/**
	 * Prints the generator matrices in standard form for dimensions 1 to @f$s@f$.
	 */
	public void printGeneratorMatrices(int s) {
		// row l, column c, dimension j.
		for (int j = 0; j < s; j++) {
			System.out.println("dim = " + (j + 1) + PrintfFormat.NEWLINE);
			for (int l = 0; l < numRows; l++) {
				for (int c = 0; c < numCols; c++)
					System.out.print(genMat[j * numCols + c][l] + "  ");
				System.out.println("");
			}
			System.out.println("----------------------------------");
		}
	}

	// Computes the digital expansion of $i$ in base $b$, and the digits
	// of the Gray code of $i$.
	// These digits are placed in arrays \texttt{bary} and \texttt{gray},
	// respectively, and the method returns the number of digits in the
	// expansion. The two arrays are assumed to be of sizes larger or
	// equal to this new number of digits.
	protected int intToDigitsGray(int b, int i, int numDigits, int[] bary, int[] gray) {
		if (i == 0)
			return 0;
		int idigits = 0; // Num. of digits in b-ary and Gray repres.
		int c;
		// convert i to b-ary representation, put digits in bary[].
		for (c = 0; i > 0; c++) {
			idigits++;
			bary[c] = i % b;
			i = i / b;
		}
		// convert b-ary representation to Gray code.
		gray[idigits - 1] = bary[idigits - 1];
		int diff;
		for (c = 0; c < idigits - 1; c++) {
			diff = bary[c] - bary[c + 1];
			if (diff < 0)
				gray[c] = diff + b;
			else
				gray[c] = diff;
		}
		for (c = idigits; c < numDigits; c++)
			gray[c] = bary[c] = 0;
		return idigits;
	}

	// ************************************************************************

	protected class DigitalNetIterator extends PointSet.DefaultPointSetIterator {

		protected int idigits;        // Num. of digits in current code.
		protected int[] bdigit;       // b-ary code of current point index.
		protected int[] gdigit;       // Gray code of current point index.
		protected int dimS;           // = dim except in the shift iterator
		                              // children where it is = dim + 1.
		protected int[] cachedCurPoint; // Digits of coords of the current point,
		// with digital shift already applied.
		// Digit l of coord j is at pos. [j*outDigits + l].
		// Has one more dimension for the shift iterators.

		public DigitalNetIterator() {
			// EpsilonHalf = 0.5 / Math.pow ((double) b, (double) outDigits);
			bdigit = new int[numCols];
			gdigit = new int[numCols];
			dimS = dim;
			cachedCurPoint = new int[(dim + 1) * outDigits];
			init();  // This call is important so that subclasses do not
			         // automatically call 'resetCurPointIndex' at the time
			         // of construction as this may cause a subtle
			         // 'NullPointerException'
		}

		public void init() { // See constructor
			resetCurPointIndex();
		}

		// We want to avoid generating 0 or 1
		public double nextDouble() {
			return nextCoordinate() + EpsilonHalf;
		}

		public double nextCoordinate() {
			if (curPointIndex >= numPoints || curCoordIndex >= dimS)
				outOfBounds();
			int start = outDigits * curCoordIndex++;
			double sum = 0;
			// Can always have up to outDigits digits, because of digital shift.
			for (int k = 0; k < outDigits; k++)
				sum += cachedCurPoint[start + k] * factor[k];
			if (digitalShift != null)
				sum += EpsilonHalf;
			return sum;
		}

		public void resetCurPointIndex() {
			if (digitalShift == null)
				for (int i = 0; i < cachedCurPoint.length; i++)
					cachedCurPoint[i] = 0;
			else {
				if (dimShift < dimS)
					addRandomShift(dimShift, dimS, shiftStream);  // Extends the random shift to
					                                              // more coordinates!
				for (int j = 0; j < dimS; j++)
					for (int k = 0; k < outDigits; k++)
						cachedCurPoint[j * outDigits + k] = digitalShift[j][k];
			}
			for (int i = 0; i < numCols; i++)
				bdigit[i] = 0;
			for (int i = 0; i < numCols; i++)
				gdigit[i] = 0;
			curPointIndex = 0;
			curCoordIndex = 0;
			idigits = 0;
		}

		public void setCurPointIndex(int i) {
			if (i == 0) {
				resetCurPointIndex();
				return;
			}
			curPointIndex = i;
			curCoordIndex = 0;
			if (digitalShift != null && dimShift < dimS)
				addRandomShift(dimShift, dimS, shiftStream);

			// Digits of Gray code, used to reconstruct cachedCurPoint.
			idigits = intToDigitsGray(b, i, numCols, bdigit, gdigit);
			int c, j, l, sum;
			for (j = 0; j < dimS; j++) {
				for (l = 0; l < outDigits; l++) {
					if (digitalShift == null)
						sum = 0;
					else
						sum = digitalShift[j][l];
					if (l < numRows)
						for (c = 0; c < idigits; c++)
							sum += genMat[j * numCols + c][l] * gdigit[c];
					cachedCurPoint[j * outDigits + l] = sum % b;
				}
			}
		}

		public int resetToNextPoint() {
			// incremental computation.
			curPointIndex++;
			curCoordIndex = 0;
			if (curPointIndex >= numPoints)
				return curPointIndex;

			// Update the digital expansion of i in base b, and find the
			// position of change in the Gray code. Set all digits == b-1 to 0
			// and increase the first one after by 1.
			int pos;      // Position of change in the Gray code.
			for (pos = 0; gdigit[pos] == b - 1; pos++)
				gdigit[pos] = 0;
			gdigit[pos]++;

			// Update the cachedCurPoint by adding the column of the gener.
			// matrix that corresponds to the Gray code digit that has changed.
			// The digital shift is already incorporated in the cached point.
			int j, l;
			int lsup = numRows;        // Max index l
			if (outDigits < numRows)
				lsup = outDigits;
			for (j = 0; j < dimS; j++) {
				for (l = 0; l < lsup; l++) {
					cachedCurPoint[j * outDigits + l] += genMat[j * numCols + pos][l];
					cachedCurPoint[j * outDigits + l] %= b;
				}
			}
			return curPointIndex;
		}
	}

	// ************************************************************************

	protected class DigitalNetIteratorNoGray extends DigitalNetIterator {

		public DigitalNetIteratorNoGray() {
			super();
		}

		public void setCurPointIndex(int i) {
			if (i == 0) {
				resetCurPointIndex();
				return;
			}
			curPointIndex = i;
			curCoordIndex = 0;
			if (dimShift < dimS)
				addRandomShift(dimShift, dimS, shiftStream);

			// Convert i to b-ary representation, put digits in bdigit.
			idigits = intToDigitsGray(b, i, numCols, bdigit, gdigit);
			int c, j, l, sum;
			for (j = 0; j < dimS; j++) {
				for (l = 0; l < outDigits; l++) {
					if (digitalShift == null)
						sum = 0;
					else
						sum = digitalShift[j][l];
					if (l < numRows)
						for (c = 0; c < idigits; c++)
							sum += genMat[j * numCols + c][l] * bdigit[c];
					cachedCurPoint[j * outDigits + l] = sum % b;
				}
			}
		}

		public int resetToNextPoint() {
			curPointIndex++;
			curCoordIndex = 0;
			if (curPointIndex >= numPoints)
				return curPointIndex;

			// Find the position of change in the digits of curPointIndex in base
			// b. Set all digits = b-1 to 0; increase the first one after by 1.
			int pos;
			for (pos = 0; bdigit[pos] == b - 1; pos++)
				bdigit[pos] = 0;
			bdigit[pos]++;

			// Update the digital expansion of curPointIndex in base b.
			// Update the cachedCurPoint by adding 1 unit at the digit pos.
			// If pos > 0, remove b-1 units in the positions < pos. Since
			// calculations are mod b, this is equivalent to adding 1 unit.
			// The digital shift is already incorporated in the cached point.
			int c, j, l;
			int lsup = numRows;        // Max index l
			if (outDigits < numRows)
				lsup = outDigits;
			for (j = 0; j < dimS; j++) {
				for (l = 0; l < lsup; l++) {
					for (c = 0; c <= pos; c++)
						cachedCurPoint[j * outDigits + l] += genMat[j * numCols + c][l];
					cachedCurPoint[j * outDigits + l] %= b;
				}
			}
			return curPointIndex;
		}
	}

}