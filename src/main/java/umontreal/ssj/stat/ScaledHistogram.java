package umontreal.ssj.stat;

/*
 * Class:        ScaledHistogram
 * Description:  A histogram whose bin heights are arbitrary real numbers
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2016  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Pierre L'Ecuyer and Mamadou Thiongane
 * @since        December 2016
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

// import java.util.logging.Logger;
// import umontreal.ssj.util.PrintfFormat;

/**
* This class provides histograms for which the bin counts (heights of rectangles)
*  are replaced by real-valued frequencies (in `double`) and can be rescaled.  
* The histogram is over a bounded interval @f$[a,b]@f$ and has a fixed number 
* of bins of equal width @f$h@f$.
*  The frequencies can be chosen (rescaled) so that the integral of the 
*  histogram is equal to a specific value.  If this value is taken as 1, 
*  i.e., if @f$h@f$ times the sum of frequencies is equal to 1, 
*  then the histogram can be seen as a density estimator for a density whose
*  support is contained in @f$[a,b]@f$.  If part of the density is outside @f$[a,b]@f$,
*  then the integral of the estimated density over @f$[a,b]@f$ can be less than 1.
*  
*  This class also implements averaged-shifted histograms (ASH) and polygonal 
*  interpolations of histograms.  See @cite tSCO15a.
* 
* <div class="SSJ-bigskip"></div>
*/

public class ScaledHistogram {
	protected int numBins; // number of bins
	protected double m_h; // width of each bin
	protected double m_a; // left boundary of first bin
	protected double m_b; // right boundary of last bin
	protected double[] height; // rescaled counters: height[j] is the height of bin j.
	protected double integral;  // Total area under the histogram, = (b-a) x sum of heights.

	private ScaledHistogram() {}

	/**
	 * Constructs a `ScaledHistogram` over the interval @f$[a,b]@f$, which is divided into 
	 * `numBins` bins of equal width. 
	 * This constructor initializes the frequency of each bin to 0.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param numBins
	 *            number of bins
	 */
	public ScaledHistogram(double a, double b, int numBins) {
		init(a, b, numBins);
	}

	/**
	 * Constructs a `ScaledHistogram` from `hist` by normalizing the bin counts so
	 * that the integral of the histogram is equal to `integral`.
	 */
	public ScaledHistogram (TallyHistogram hist, double integral) {
		init (hist, integral);
	}


	/**
	 * Initializes the `ScaledHistogram` so it covers the interval @f$[a,b]@f$, 
	 * which is divided into `numBins` bins of equal width. 
	 * Is initializes the frequency of each bin to 0.
	 * @param a
	 *            left boundary of interval
	 * @param b
	 *            right boundary of interval
	 * @param numBins
	 *            number of bins
	 */
	public void init (double a, double b, int numBins) {
		if (b <= a)
			throw new IllegalArgumentException("   b <= a");
		this.numBins = numBins;
		m_h = (b - a) / numBins;
		m_a = a;
		m_b = b;
		height = new double[numBins];
		for (int i = 0; i < numBins; i++)
			height[i] = 0;
		integral = 0.0;
	}

	/**
	 * Initializes this `ScaledHistogram` using the @ref TallyHistogram `hist`.
	 * It uses the same interval @f$[a,b]@f$, same bins, and rescaled the counters
	 * so the integral of the histogram equal the value specified by `integral`.
	 */
	public void init (TallyHistogram hist, double integral) {
		m_a = hist.getA();
		m_b = hist.getB();
		m_h = hist.getH();
		numBins = hist.numBins;
		height = new double[numBins];
		this.integral = integral;
		int count[] = hist.getCounters();
		double scaleFactor = integral / (hist.numberObs() * m_h);
		for (int i = 0; i < numBins; i++) 
			height[i] = count[i] * scaleFactor;  
	}
	
	/**
	 * Initializes all the heights (frequencies) to 0.
	 */
	public void init () {
		for (int i = 0; i < numBins; i++)
			height[i] = 0.0;
	}
		
	/**
	 * Rescales the histogram by renormalizing the frequencies so its integral has the value
	 * specified by `integral`.
	 */
	public void rescale (double integral) {
		double scaleFactor = integral / this.integral;
		for (int i = 0; i < numBins; i++)
			height[i] *= scaleFactor;
		this.integral = integral;
	}

	/**
	 * Returns an ASH-transformed version of this scaled histogram. The
	 * ASH-transformed histogram has the same bin size as the original. The new
	 * frequency (height) in any given bin is the weighted average of the frequencies in
	 * the neighboring bins, with weights @f$(r-d)/r^2@f$ given to bins that are at
	 * distance @f$d@f$ from the target bin, for all @f$d < r@f$.
	 * 
	 */
	public ScaledHistogram averageShiftedHistogram (int r) {
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
		for (int k = 0; k < numBins; k++) {
			heightNew[k] = r * height[k]; 
			for (int ell = 1; ell < r; ell++) {
				if (k-ell >= 0) heightNew[k] += (r-ell) * height[k-ell];
				if (k+ell < numBins) heightNew[k] += (r-ell) * height[k+ell];	
			}		
			heightNew[k] *= rscale;
			sum += heightNew[k];
   	    }
		image.height = heightNew;
		image.integral = sum*m_h;
		return image;
	}

	/**
	 * Similar to `averageShiftedHistogram`, except that it assumes that the density 
	 * is over a close interval @f$[a,b]@f$ and is rescaled differently for the intervals 
	 * that are near the boundary, to account for the fact that the intervals 
	 * outside the boundaries are not counted.
	 */
	public ScaledHistogram averageShiftedHistogramTrunc (int r) {
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		// double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
		for (int k = 0; k < numBins; k++) {
			heightNew[k] = r * height[k];  // 
			int weight = r;
			for (int ell = 1; ell < r; ell++) {
				if (k-ell >= 0) {
					heightNew[k] += (r-ell) * height[k-ell];
					weight += r-ell;
				}
				if (k+ell < numBins) {
					heightNew[k] += (r-ell) * height[k+ell];	
					weight += r-ell;
				}
			}		
			heightNew[k] *= 1.0 / (double)weight;
			sum += heightNew[k];
   	    }
		image.height = heightNew;
		image.integral = sum * m_h;
		return image;
	}

	/**
	 * Similar to `averageShiftedHistogram`, except that uses a weighted average.
	 * For the new average in a given bin, any neighbor bin at distance @f$\ell < r@f$ is given 
	 * a weight proportional to `w[i]`. The given weights do not have to sum to 1;
	 * they are rescaled so the sum of weights that go to any given bin is 1.
	 */
	public ScaledHistogram averageShiftedHistogram (int r, double[] w) {
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		double weight = w[0];
		for (int ell = 1; ell < r; ell++)
            weight += 2.0 * w[ell];
	    double rscale = 1.0 / weight;     // Rescaling factor for each bin.
		double sum = 0.0;
		for (int k = 0; k < numBins; k++) {
			heightNew[k] = w[0] * height[k];  // 
			for (int ell = 1; ell < r; ell++) {
				if (k-ell >= 0) {
					heightNew[k] += w[ell] * height[k-ell];
				}
				if (k+ell < numBins) {
					heightNew[k] += w[ell] * height[k+ell];	
				}
			}		
			heightNew[k] *= rscale;
			sum += heightNew[k];
   	    }
		image.height = heightNew;
		image.integral = sum * m_h;
		return image;
	}
		
		
	/**
	 * Similar to `averageShiftedHistogramTrunc`, except that uses a weighted average.
	 * For the new average in a given bin, any neighbor bin at distance @f$\ell < r@f$ is given 
	 * a weight proportional to `w[i]`. The given weights do not have to sum to 1;
	 * they are rescaled so the sum of weights that go to any given bin is 1 (not counting 
	 * the weights given to bins that fall outside the interval).
	 */
	public ScaledHistogram averageShiftedHistogramTrunc (int r, double[] w) {
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		// double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
		for (int k = 0; k < numBins; k++) {
			heightNew[k] = w[0] * height[k];  // 
			double weight = w[0];
			for (int ell = 1; ell < r; ell++) {
				if (k-ell >= 0) {
					heightNew[k] += w[ell] * height[k-ell];
					weight += w[ell];
				}
				if (k+ell < numBins) {
					heightNew[k] += w[ell] * height[k+ell];	
					weight += w[ell];
				}
			}		
			heightNew[k] *= 1.0 / weight;
			sum += heightNew[k];
   	    }
		image.height = heightNew;
		image.integral = sum * m_h;
		return image;
	}

	/*	public ScaledHistogram averageShiftedHistogram1 (int r) {
		// if (numBins % r != 0)
			// throw new IllegalArgumentException("r must divide the number of bins.");
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
        double S1[] = new double[numBins];
        double S2[] = new double[numBins];
        S1[0] = height[0];
        S2[0] = 0.0;
        heightNew[0] = r * S1[0];
		for (int ell = 1; ell < Math.min(r, numBins); ell++) {
			S2[0] += height[ell];
			heightNew[0] += (r-ell) * height[ell];
		}
       	for (int k = 1; k < numBins; k++) {
      		for (int ell = 1; ell < r; ell++) {
      			S1[k] = S1[k-1] + height[k];
      			if (k >= r) S1[k] -= height[k-r];
      			S2[k] = S2[k-1] - height[k];
      			if (k+r-1 < numBins) S2[k] += height[k+r-1];	
      			heightNew[k] = heightNew[k-1] + S2[k-1] - S1[k-1];
			}
			heightNew[k] *= rscale;
			sum += heightNew[k];
   	    }
		image.height = heightNew;
		image.integral = sum*m_h;
		return image;
	}*/
	

	/**
	 * This is supposed to be a faster implementation of `averageShiftedHistogram(r)`.
	 */
	public ScaledHistogram averageShiftedHistogram1 (int r) {
		ScaledHistogram image = clone();
		double[] heightNew = image.getHeights();
		double rscale = 1.0 / (r * r);   // Rescaling to be made for each bin.
		double sum = 0.0;
        double S1[] = new double[numBins];
        double S2[] = new double[numBins];
        S1[0] = height[0];
        S2[0] = 0.0;
        heightNew[0] = r * S1[0];
		for (int ell = 1; ell <= Math.min(r, numBins); ell++){
			S2[0] += height[ell];	
			heightNew[0] += (r-ell) * height[ell];
		}
	    for (int k = 2; k <= numBins; k++) {
      			S1[k-1] = S1[k-2] + height[k-1];
      			if (k >= r) S1[k-1] -= height[k-r];
      			S2[k-1] = S2[k-2] - height[k-1];
      			if (k+r< numBins) S2[k-1] += height[k+r-1];	
      			heightNew[k-1] = heightNew[k-2] + S2[k-2] - S1[k-2];
		}      	
    	for (int k = 0; k <numBins; k++) {
    		heightNew[k] *= rscale;
			sum += heightNew[k];	
    	}
		image.height = heightNew;
		image.integral = sum*m_h;
		return image;
	}
	
  
     /**
	 * Returns the number of bins @f$s@f$ dividing the interval 
	 * @f$[a,b]@f$. Does not count the two extra bins for the values of
	 * @f$x<a@f$ or @f$x>b@f$.
	 * @return the number of bins
	 */
	public int getNumBins() {
		return numBins;
	}

	/**
	 * Returns the left boundary @f$a@f$ of interval @f$[a,b]@f$.
	 * 
	 * @return left boundary of interval
	 */
	public double getA() {
		return m_a;
	}

	/**
	 * Returns the right boundary @f$b@f$ of interval @f$[a,b]@f$.
	 * 
	 * @return right boundary of interval
	 */
	public double getB() {
		return m_b;
	}

	/**
	 * return the array counts of the histogram.
	 */
	public double[] getHeights() {
		return height;
	}

	/**
	 * return the integral the histogram.
	 */
	public double getIntegral() {
		return integral;
	}

	/**
	 * Computes and returns the integrated square error (ISE) of a histogram 
	 * w.r.t. the @f$U(0,1)@f$ distribution.  Assumes the histogram integrates to 1.
	 */
	public double ISEvsU01 () {
		double sum = 0.0;
		for (int j = 0; j < numBins; ++j)
		   sum += (height[j] - 1.0) * (height[j] - 1.0);
		return sum / numBins;
	}
	
	/**
	 * Computes and returns the ISE of a polygonal density w.r.t. the @f$U(0,1)@f$
	 * distribution.   Assumes the histogram integrates to 1.
	 * Over interval @f$(a,b)@f$ = '(w0[j],w0[j+1])', the ISE is 
	 * @f$\int_0^1 ((a + (b-a)x)^2 dx = (b^3 - a^3) / (3 (b-a)) = (b^2 + a^2 + ab)/3.@f$
	 */
	public double ISEvsU01polygonal () {
		double w0[] = new double[numBins];  // histogram shifted down to zero mean.
		for (int j = 0; j < numBins; ++j) {
			w0[j] = height[j] - 1.0;
		}
		double a, b;
		double sum = 0.5 * (w0[0] * w0[0] + w0[numBins-1] * w0[numBins-1]);
		for (int j = 0; j < numBins-2; ++j) {
		   a = w0[j];
		   b = w0[j+1];
		   sum += 0.3333333333333333333 * (b*b + a*a + a*b);
		}
		return sum / (double)numBins;
	}
	
	/**
	 * Clones this object and the array which stores the counters.
	 */
	public ScaledHistogram clone() {
		ScaledHistogram image = new ScaledHistogram();
		image.numBins = numBins;
		image.m_h = m_h;
		image.m_a = m_a;
		image.m_b = m_b;
		image.height = new double[numBins]; 
		image.integral = integral;  
		for (int j = 1; j < numBins; ++j)
		   image.height[j] = height[j];
		return image;
	}
		
}
