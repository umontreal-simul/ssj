package umontreal.ssj.stat;

/*
 * Class:        HistogramOnly
 * Description:  Histogram of a tally
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Pierre L'Ecuyer
 * @since        January 2017
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
// package umontreal.ssj.stat;
import umontreal.ssj.stat.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class is similar to @ref TallyHistogram, except that it does not maintain the min, max,
 * average, and variance of the observations. Only the counters for the histogram are maintained.
 * Also it does not maintain counters for virtual bins on the left and on the right to count 
 * the observations that fall outside @f$[a,b]@f$.
 * The methods that are supposed to return these values return an error message instead.
 * The only advantage of not maintaining these counters and values is to increase speed (slightly). 
 * We define this class as an extension of @ref TallyHistogram to avoid duplicating code.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class HistogramOnly extends TallyHistogram implements Cloneable {

	/**
	 * Constructs a `HistogramOnly` statistical probe. Divide the interval
	 * 
	 * @f$[a,b]@f$ into @f$s@f$ bins of equal width and initializes a counter to 0 for each bin.
	 *             Whenever an observation falls into a bin, the bin counter is increased by 1.
	 *             There are two extra bins (and counters) that count the number of
	 *             observations @f$x@f$ that fall outside the interval @f$[a,b]@f$: one for
	 *             those @f$x< a@f$, and the other for those @f$x > b@f$.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param s
	 *            number of bins
	 */
	public HistogramOnly(double a, double b, int s) {
		super(a, b, s);
	}

	/**
	 * Constructs a new `HistogramOnly` statistical probe with name `name`.
	 * 
	 * @param name
	 *            the name of the tally.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param s
	 *            number of bins
	 */
	public HistogramOnly(String name, double a, double b, int s) {
		super(a, b, s);
		this.name = name;
	}


	@Override
	public void init() {
		numObs = 0;
		for (int i = 0; i < numBins; i++)
			count[i] = 0;
	}
	
	/**
	 * Gives a new observation @f$x@f$ to the statistical collectors. Increases by 1 the bin counter
	 * in which value @f$x@f$ falls. 
	 * 
	 * @param x
	 *            observation value
	 */
	public void add(double x) {
		numObs++;
		if ((x >= m_a) & (x <= m_b))
			++count[(int) ((x - m_a) / m_h)];
	}

	   @Override
	   public double sum() {
			throw new IllegalStateException("HistogramOnly.variance() is not supported.");
	 	    // return Double.NaN;
	   }

	   @Override
	   public double average() {
			throw new IllegalStateException("HistogramOnly.variance() is not supported.");
	 	    // return Double.NaN;
	   }

	   @Override
	   public double variance() {
			throw new IllegalStateException("HistogramOnly.variance() is not supported.");
	 	    // return Double.NaN;
	   }


	/**
	 * Clones this object and the array which stores the counters.
	 */
	public HistogramOnly clone() {
		HistogramOnly image = (HistogramOnly) super.clone();
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
		for (int i = 0; i < numBins; i++) {
			double a = m_a + (i - 1) * m_h;
			double b = m_a + i * m_h;
			sb.append("   (" + PrintfFormat.f(6, 3, a) + ", " + PrintfFormat.f(6, 3, b) + ")    "
			        + count[i] + PrintfFormat.NEWLINE);
		}
		sb.append("}" + PrintfFormat.NEWLINE);
		return sb.toString();
	}

	@Override
	public String report() {
		// TODO Auto-generated method stub
		throw new IllegalStateException("HistogramOnly.report() is not supported.");
	}

	@Override
	public String shortReport() {
		// TODO Auto-generated method stub
		throw new IllegalStateException("HistogramOnly.shortReport() is not supported.");
	}

	@Override
	public String shortReportHeader() {
		// TODO Auto-generated method stub
		throw new IllegalStateException("HistogramOnly.shortReportHeader() is not supported.");
	}

}