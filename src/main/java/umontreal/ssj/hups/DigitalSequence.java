/*
 * Class:        DigitalSequence
 * Description:  abstract class with methods specific to digital sequences
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

/**
 * This abstract class describes methods specific to digital sequences. Concrete classes must
 * implement the #extendSequence method that increases the number of points of the digital sequence.
 * Calling the methods #toNet or #toNetShiftCj will transform the digital sequence into a digital
 * net, which has a fixed number of points @f$n@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class DigitalSequence extends DigitalNet {

	/**
	 * Increases the number of points to @f$n = b^k@f$ from now on.
	 * 
	 * @param k
	 *            there will be b^k points
	 */
	public abstract void extendSequence(int k);

	private int[][] copyDigitalShift(int[][] S) {
		// Copy the shift S into T and returns T.
		if (S == null)
			return null;
		int[][] T = new int[S.length][S[0].length];
		for (int i = 0; i < S.length; ++i)
			for (int j = 0; j < S[0].length; ++j)
				T[i][j] = S[i][j];
		return T;
	}

	private DigitalNet initNetVar(boolean shiftFlag) {
		// Initializes the net for the two toNet methods below.
		DigitalNet net = new DigitalNet();
		if (shiftFlag)
			net.dim = dim + 1;
		else
			net.dim = dim;
		net.numPoints = numPoints;
		net.numCols = numCols;
		net.numRows = numRows;
		net.outDigits = outDigits;
		net.normFactor = normFactor;
		net.b = b;
		net.factor = new double[outDigits];
		for (int i = 0; i < outDigits; i++)
			net.factor[i] = factor[i];
		net.genMat = new int[net.dim * numCols][numRows];
		net.shiftStream = shiftStream;
		net.capacityShift = capacityShift;
		net.dimShift = dimShift;
		net.digitalShift = copyDigitalShift(digitalShift);
		if (shiftFlag && shiftStream != null) {
			net.addRandomShift(dimShift, dimShift + 1, shiftStream);
		}
		return net;
	}

	/**
	 * Transforms this digital sequence into a digital net without changing the coordinates of the
	 * points. Returns the digital net.
	 */
	public DigitalNet toNet() {
		DigitalNet net = initNetVar(false);
		final int N = dim * numCols;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < numRows; j++)
				net.genMat[i][j] = genMat[i][j];
		return net;
	}

	/**
	 * Transforms this digital sequence into a digital net by adding one dimension and shifting all
	 * coordinates by one position. The first coordinate of point @f$i@f$ is @f$i/n@f$,
	 * where @f$n@f$ is the total number of points. Thus if the coordinates of a point of the
	 * digital sequence were @f$(x_0, x_1, x_2, …, x_{s-1})@f$, then the coordinates of the point of
	 * the digital net will be @f$(i/n, x_0, x_1, …, x_{s-1})@f$. In other words, for the digital
	 * net,
	 * 
	 * @f$\mathbf{C}_0@f$ is the reflected identity and for @f$j\ge1@f$, the @f$\mathbf{C}_j@f$ used
	 *                    is the @f$\mathbf{C}_{j-1}@f$ of the digital sequence. If the digital
	 *                    sequence uses a digital shift, then the digital net will include the
	 *                    digital shift with one more dimension also. Returns the digital net.
	 */
	public DigitalNet toNetShiftCj() {
		DigitalNet net = initNetVar(true);
		int j, c, l, start;

		/* Shift all coordinates from the sequence by 1 dimension */
		for (j = dim; j >= 1; j--) {
			start = j * numCols;
			for (c = 0; c < numCols; c++)
				for (l = 0; l < numRows; l++)
					net.genMat[start + c][l] = genMat[start - numCols + c][l];
		}

		// j = 0: initialize C_0 to the reflected identity.
		for (c = 0; c < numCols; c++) {
			for (l = 0; l < numRows; l++)
				net.genMat[c][l] = 0;
			net.genMat[c][numCols - c - 1] = 1;
		}
		return net;
	}

	/**
	 * Similar to #iterator, except that the first coordinate of the points is @f$i/n@f$, the second
	 * coordinate is obtained via the generating matrix @f$\mathbf{C}_0@f$, the next one via
	 * 
	 * @f$\mathbf{C}_1@f$, and so on. Thus, this iterator shifts all coordinates of each point one
	 *                     position to the right and sets the first coordinate of point @f$i@f$
	 *                     to @f$i/n@f$, so that the points enumerated with this iterator have one
	 *                     more dimension. A digital shift, if present, will have one more dimension
	 *                     also. This iterator uses the Gray code.
	 */
	public PointSetIterator iteratorShift() {
		return new DigitalNetIteratorShiftGenerators();
	}

	/**
	 * This iterator shifts all coordinates of each point one position to the right and sets the
	 * first coordinate of point @f$i@f$ to
	 * 
	 * @f$i/n@f$, so that the points enumerated with this iterator have one more dimension. This
	 *            iterator does not use the Gray code; the points are enumerated in the order of
	 *            their first coordinate before randomization. A digital shift, if present, will
	 *            have one more dimension also.
	 */
	public PointSetIterator iteratorShiftNoGray() {
		return new DigitalNetIteratorShiftNoGray();
	}

	// ************************************************************************

	protected class DigitalNetIteratorShiftGenerators extends DigitalNetIterator {
		// Similar to DigitalNetIterator; the first coordinate
		// of point i is i/n, and all the others are shifted one position
		// to the right. The points have dimension = dim + 1.

		public DigitalNetIteratorShiftGenerators() {
			super();
			dimS = dim + 1;
			if (digitalShift != null && dimShift < dimS)
				addRandomShift(dimShift, dimS, shiftStream);
			init2();
		}

		public void init() {   // This method is necessary to overload
		}                      // the init() of DigitalNetIterator

		public void init2() { // See constructor
			resetCurPointIndex();
		}

		public void setCurPointIndex(int i) {
			if (i == 0) {
				resetCurPointIndex();
				return;
			}
			curPointIndex = i;
			curCoordIndex = 0;

			// Digits of Gray code, used to reconstruct cachedCurPoint.
			idigits = intToDigitsGray(b, i, numCols, bdigit, gdigit);
			int c, j, l, sum;
			for (j = 1; j <= dim; j++) {
				for (l = 0; l < outDigits; l++) {
					if (digitalShift == null)
						sum = 0;
					else
						sum = digitalShift[j][l];
					if (l < numRows)
						for (c = 0; c < idigits; c++)
							sum += genMat[(j - 1) * numCols + c][l] * gdigit[c];
					cachedCurPoint[j * outDigits + l] = sum % b;
				}
			}
			// The case j = 0
			for (l = 0; l < outDigits; l++) {
				if (digitalShift == null)
					sum = 0;
				else
					sum = digitalShift[0][l];
				if (l < numRows)
					for (c = 0; c < idigits; c++)
						if (l == numCols - c - 1)
							sum += gdigit[c];
				cachedCurPoint[l] = sum % b;
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
			int c, j, l;
			int lsup = numRows;        // Max index l
			if (outDigits < numRows)
				lsup = outDigits;
			for (j = 1; j <= dim; j++) {
				for (l = 0; l < lsup; l++) {
					cachedCurPoint[j * outDigits + l] += genMat[(j - 1) * numCols + pos][l];
					cachedCurPoint[j * outDigits + l] %= b;
				}
			}
			// The case j = 0
			l = numCols - pos - 1;
			if (l < lsup) {
				cachedCurPoint[l] += 1;
				cachedCurPoint[l] %= b;
			}

			return curPointIndex;
		}
	}

	// ************************************************************************

	protected class DigitalNetIteratorShiftNoGray extends DigitalNetIterator {
		// Similar to DigitalNetIterator; the first coordinate
		// of point i is i/n, and all the others are shifted one position
		// to the right. The points have dimension = dim + 1.

		public DigitalNetIteratorShiftNoGray() {
			super();
			dimS = dim + 1;
			if (digitalShift != null && dimShift < dimS)
				addRandomShift(dimShift, dimS, shiftStream);
			init2();
		}

		public void init() {   // This method is necessary to overload
		}                      // the init() of DigitalNetIterator

		public void init2() { // See constructor
			resetCurPointIndex();
		}

		public void setCurPointIndex(int i) {
			if (i == 0) {
				resetCurPointIndex();
				return;
			}
			curPointIndex = i;
			curCoordIndex = 0;

			// Convert i to b-ary representation, put digits in bdigit.
			idigits = intToDigitsGray(b, i, numCols, bdigit, gdigit);
			int c, j, l, sum;
			for (j = 1; j <= dim; j++) {
				for (l = 0; l < outDigits; l++) {
					if (digitalShift == null)
						sum = 0;
					else
						sum = digitalShift[j][l];
					if (l < numRows)
						for (c = 0; c < idigits; c++) {
							sum += genMat[(j - 1) * numCols + c][l] * bdigit[c];
							sum %= b;
						}
					cachedCurPoint[j * outDigits + l] = sum;
				}
			}
			// The case j = 0
			for (l = 0; l < outDigits; l++) {
				if (digitalShift == null)
					sum = 0;
				else
					sum = digitalShift[0][l];
				if (l < numRows)
					for (c = 0; c < idigits; c++)
						if (l == numCols - c - 1)
							sum += bdigit[c];
				cachedCurPoint[l] = sum % b;
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
			for (j = 1; j <= dim; j++) {
				for (l = 0; l < lsup; l++) {
					for (c = 0; c <= pos; c++)
						cachedCurPoint[j * outDigits + l] += genMat[(j - 1) * numCols + c][l];
					cachedCurPoint[j * outDigits + l] %= b;
				}
			}
			// The case j = 0
			for (l = 0; l < lsup; l++) {
				for (c = 0; c <= pos; c++)
					if (l == numCols - c - 1) {
						cachedCurPoint[l] += 1;
						cachedCurPoint[l] %= b;
					}
			}

			return curPointIndex;
		}

	}
}