package umontreal.ssj.stat;

/*
 * Class:        TallyHistogram
 * Description:  A Tally that also builds a Histogram
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard and Pierre L'Ecuyer
 * @since        January 2011
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

import umontreal.ssj.util.PrintfFormat;

/**
 * This class extends @ref Tally. It does not store individual observations, but in addition to 
 * maintaining the counters as in @ref Tally, it also constructs a histogram for the observations. 
 * The histogram is over a bounded interval @f$[a,b]@f$ and has a fixed number of bins of equal width,
 * both specified by the user. The number of observations falling into each of the bins 
 * is kept in an array of counters. This array can be accessed directly by the user.
 * Note that one should never add or remove observations *directly* on this array of 
 * bin counters because this would put the @ref Tally counters in an inconsistent state.
 * Additional variables count the number of observations falling outside the interval @f$[a,b]@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class TallyHistogram extends Tally {
	protected int numBins; // number of bins
	protected int[] count; // counter for number of values in bin[i].
	protected int leftCount;  // Count values that are less than a.
	protected int rightCount; // Count values that are larger than b.
	protected double m_h; // width of 1 bin
	protected double m_a; // left boundary of first bin
	protected double m_b; // right boundary of last bin

	/**
	 * Constructs a `TallyHistogram` statistical probe. Divide the interval
	 * 
	 * @f$[a,b]@f$ into @f$s@f$ bins of equal width and initializes a counter to 0 for each bin.
	 *             Whenever an observation falls into a bin, the bin counter is increased by 1.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param numBins
	 *            number of bins (of equal width)
	 */
	public TallyHistogram(double a, double b, int numBins) {
		super();
		init(a, b, numBins);
	}

	/**
	 * Constructs a new `TallyHistogram` statistical probe with name `name`.
	 * 
	 * @param name
	 *            the name of the tally.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param numBins
	 *            number of bins
	 */
	public TallyHistogram(String name, double a, double b, int numBins) {
		super(name);
		init(a, b, numBins);
	}

	/**
	 * Initializes this object. Divide the interval @f$[a,b]@f$ into
	 * 
	 * @f$s@f$ bins of equal width and initializes all counters to 0.
	 * @param numBins
	 *            number of bins
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 */
	public void init(double a, double b, int numBins) {
		/*
		 * The counters count[0] to count[s-1] contains the number of observations falling in the
		 * interval [a, b]. leftCount is the number of observations < a, and rightCount is the
		 * number of observations > b.
		 */
		super.init();
		if (b <= a)
			throw new IllegalArgumentException("   b <= a");
		count = new int[numBins];
		this.numBins = numBins;
		m_h = (b - a) / (double)numBins;
		m_a = a;
		m_b = b;
		leftCount = rightCount = 0;
		for (int i = 0; i < numBins; i++)
			count[i] = 0;
	}

	/**
	 * Initializes all the counters and accumulators, including those of the `Tally` object.
	 *
	 */
	public void init() {
		super.init();
		leftCount = rightCount = 0;
		for (int i = 0; i < numBins; i++)
			count[i] = 0;
	}

	
	/**
	 * Fills this object from the first numObs observations in array obs.
	 */
	public void fillFromArray(double[] obs, int numObs) {
		init();
		for (int i = 0; i < numObs; i++)
			add(obs[i]);
	}

	/**
	 * Fills this object from the entire array obs.
	 */
	public void fillFromArray(double[] obs) {
		fillFromArray(obs, obs.length);
	}

	/**
	 * Fills this object from the observations in a TallyStore object.
	 */
	public void fillFromTallyStore(TallyStore ts) {
		fillFromArray(ts.getArray(), ts.numberObs());
	}

	/**
	 * Gives a new observation @f$x@f$ to the statistical probe. Updates are made as for the
	 * parent `Tally` object.  Also increases by 1 the bin counter
	 * in which value @f$x@f$ falls. Values that fall outside the interval @f$[a,b]@f$ are added 
	 * to the extra bin counters.
	 * 
	 * @param x
	 *            observation value
	 */
	public void add(double x) {
		super.add(x);
		if (x < m_a)
			++leftCount;
		else if (x > m_b)
			++rightCount;
		else {
			int i = (int) ((x - m_a) / m_h);
			++count[i];
		}
	}

	/**
	 * Remove empty bins in the tails (left and right), without changing the bin size.
	 * This gives a new @ref TallyHistogram which may have fewer bins.
	 */
	public TallyHistogram trimHistogram() {
		TallyHistogram image = (TallyHistogram) super.clone();
		int i = 0;
		int j = numBins-1; // last bin in the initial histogram
		int cpL = 0; // number of empty bins from left initialized to zero
		int cpR = 0; // number of empty bins from right initialized to zero
		while (count[i] == 0) {
			i++;
			cpL++;
		}
		while (count[j] == 0) {
			j--;
			cpR++;
		}
		int[] coco = new int[numBins - cpL - cpR];
		System.arraycopy(count, i, coco, 0, j - i + 1);
		image.count = coco;
		image.m_h = m_h;
		image.m_a = m_a + (cpL * m_h);
		image.m_b = m_b - (cpR * m_h);
		image.numBins = numBins - cpL - cpR;
		image.leftCount = leftCount;
		image.rightCount = rightCount;
		return image;
	}

	/**
	 * Merges this histogram with the other histogram, by adding the bin counts of the two
	 * histograms.
	 * 
	 * @param other
	 *            the histogram to add
	 * 
	 *            Returns the merged histogram.
	 */
	public TallyHistogram addHistograms(TallyHistogram other) {
		if (this.numBins != other.numBins)
			throw new IllegalArgumentException("different number of bin in two histogram to merge");
		TallyHistogram image = (TallyHistogram) super.clone();
		int[] countNew = new int[numBins];
		System.arraycopy(count, 0, countNew, 0, numBins);
		int coOther[] = other.getCounters();
		for (int i = 0; i < countNew.length; i++)
			countNew[i] = countNew[i] + coOther[i];
		image.count = countNew;
		image.leftCount = leftCount + other.leftCount;
		image.rightCount = rightCount + other.rightCount;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.numBins = numBins;
		return image;

	}

	/**
	 * Merges bins by groups of size @f$g@f$. If there are @f$m@f$ bins initially, the new number of bins
	 * will be @f$\lceil m/g\rceil@f$. The last bin may regroup less than @f$g@f$ original bins 
	 * if @f$m@f$ is not a multiple of @f$g@f$. In this case the upper bound @f$b@f$ is increased accordingly.
	 **/
	public TallyHistogram aggregateBins(int g) {
		TallyHistogram image = (TallyHistogram) super.clone();
		int numBinsNew = (int)Math.ceil((double) numBins / (double) g);
		int[] countNew = new int[numBinsNew];
		int b = 0;
		for (int j = 0; j < numBinsNew-1; j++) {
			for (int i = b; i < b + g; i++)
				countNew[j] += count[i];
			b = b + g;
		}
		while (b < numBins - 1) {
			countNew[numBinsNew-1] += count[b];
			b++;
		}
		image.count = countNew;
		image.m_h = m_h * g;
		image.m_a = m_a;
		image.m_b = m_h * numBinsNew;
		image.numBins = numBinsNew;
		image.leftCount = leftCount;
		image.rightCount = rightCount;
		return image;
	}

	/**
	 * Returns the array of bin counters. Each counter contains the number of observations that fell in its
	 * corresponding bin. 
	 * 
	 * @return the array of bin counters
	 */
	public int[] getCounters() {
		return count;
	}

	/**
	 * Returns the number of bins @f$s@f$.  
	 * 
	 * @return the number of bins
	 */
	public int getNumBins() {
		return numBins;
	}

	/**
	 * Returns the left boundary @f$a@f$ of the interval @f$[a,b]@f$.
	 * 
	 * @return left boundary of interval
	 */
	public double getA() {
		return m_a;
	}

	/**
	 * Returns the right boundary @f$b@f$ of the interval @f$[a,b]@f$.
	 * 
	 * @return right boundary of interval
	 */
	public double getB() {
		return m_b;
	}

	/**
	 * Returns the width @f$h@f$ of the bins.
	 * 
	 * @return the width of the bins
	 */
	public double getH() {
		return m_h;
	}

	/**
	 * Returns the proportion of the collected observations that lie within the boundaries 
	 * @f$[a,b]@f$ of the histogram;  that is, the number that fell within @f$[a,b]@f$
	 * divided by the total number that were collected.
	 * 
	 * @return the proportion of observations that lie in the histogram.
	 */
	public double getProportionInBoundaries() {
		int total = 0;
		for(int num : count) {
			total += num;
		}
		return  (double) total / (double)(total + leftCount + rightCount); 
	}
	
	/**
	 * Clones this object and the array that stores the counters.
	 */
	public TallyHistogram clone() {
		TallyHistogram image = (TallyHistogram) super.clone();
		int[] coco = new int[numBins];
		System.arraycopy(count, 0, coco, 0, numBins);
		image.count = coco;
        image.leftCount = leftCount;
        image.rightCount = rightCount;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.numBins = numBins;
		return image;
	}

	/**
	 * Returns the bin counters as a `String`.
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("---------------------------------------" + PrintfFormat.NEWLINE);
		sb.append(name + PrintfFormat.NEWLINE);
		sb.append("Interval = [ " + m_a + ", " + m_b + " ]" + PrintfFormat.NEWLINE);
		sb.append("Number of bins = " + numBins + PrintfFormat.NEWLINE);
		sb.append(PrintfFormat.NEWLINE + "Counters = {" + PrintfFormat.NEWLINE);
		sb.append("   (-inf, " + PrintfFormat.f(6, 3, m_a) + ")    " + leftCount
		        + PrintfFormat.NEWLINE);
		for (int i = 0; i < numBins; i++) {
			double a = m_a + (i - 1) * m_h;
			double b = m_a + i * m_h;
			sb.append("   (" + PrintfFormat.f(6, 3, a) + ", " + PrintfFormat.f(6, 3, b) + ")    "
			        + count[i] + PrintfFormat.NEWLINE);
		}
		sb.append("   (" + PrintfFormat.f(6, 3, m_b) + ", inf)    " + rightCount
		        + PrintfFormat.NEWLINE);
		sb.append("}" + PrintfFormat.NEWLINE);
		return sb.toString();
	}

}