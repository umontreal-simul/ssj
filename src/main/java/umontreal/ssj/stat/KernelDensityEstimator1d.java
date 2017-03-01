package umontreal.ssj.stat;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
/**
Uses Apache Commons Math Library (tested with Apache Commons Math 3.6.1)
http://commons.apache.org/proper/commons-math/

Reliable and extremely fast kernel density estimator for one-dimensional data;
Gaussian kernel is assumed and the bandwidth is chosen automatically;
Unlike many other implementations, this one is immune to problems
caused by multimodal densities with widely separated modes. The
estimation does not deteriorate for multimodal densities, because we never assume
a parametric model for the data.

Call the method "kde" to create a kernel density estimate and cdf estimate from input data.
The estimate will be stored in the output variables described below.
They can be read through accessor methods.

INPUTS:
data    - an array of data from which the density estimate is constructed;
n  - the number of mesh points used in the uniform discretization of the
interval [min, max]; n has to be a power of two; if n is not a power of two, then
n is rounded up to the next power of two, i.e., n is set to n=2^ceil(log2(n));
the default value of n is n=2^14;
min, max  - defines the interval [min,max] on which the density estimate is constructed;
the default values of min and max are:
min=min(data)-range/10 and max=max(data)+range/10, where range=max(data)-min(data);

OUTPUTS:
bandwidth - the optimal bandwidth (Gaussian kernel assumed);
density - array of length 'n' with the values of the density
estimate at the grid points;
xmesh   - the grid over which the density estimate is computed;
cdf  - array of length 'n' with the values of the cdf;
bandwidthCdf - the optimal bandwidth for cdf estimation

References:
Kernel density estimation via diffusion
Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
D. P. Kroese, T. Taimre, Z. I. Botev (2011)
Handbook of Monte Carlo Methods, pp. 319-330.

Example:
double[] data = some array of data;
KernelDensityEstimator1d kde1d = new KernelDensityEstimator1d();
kde1d.kde(data, (int)Math.pow(2, 14), -5, 5);
System.out.println(kde1d.getBandwidth());
*/
public class KernelDensityEstimator1d {
	
	private double bandwidth;
	private double[] density;
	private double[] xmesh;
	private double[] cdf;
	private double bandwidthCdf;
	
	//accessor methods to read the output variables
	public double getBandwidth() {
		return bandwidth;
	}
	
	public double[] getDensity() {
		return density;
	}
	
	public double[] getXmesh() {
		return xmesh;
	}
	
	public double[] getCdf() {
		return cdf;
	}
	
	public double getBandwidthCdf() {
		return bandwidthCdf;
	}
	
	public void kde(double[] data) {
		//if n is not supplied switch to the default
		kde(data, (int)Math.pow(2, 14));
	}
	
	public void kde(double[] data, int n) {
		//define the default interval [min,max]
		double min = data[0];
		double max = data[0];
		for (double d : data) {
			if (d<min) {
				min = d;
			}
			if (d>max) {
				max = d;
			}
		}
		double range = max-min;
		kde(data, n, min-range/10, max+range/10);
	}
	
	// The bandwidth is for the interval rescaled to [0,1].  Density is computed at n points.
	public void kde(double[] data, int n, double min, double max, double bandwidth) {
		//round up n to the next power of 2
		n = (int)Math.pow(2, Math.ceil(Math.log(n)/Math.log(2)));
		
		//set up the grid over which the density estimate is computed
		double range = max-min;
		bandwidth *= range;
		double dx = range/n;
		xmesh = new double[n];   // Points at which the density is computed.
		//  System.out.println("kde, n = " + n);

		xmesh[0] = min+dx*0.5;
		for (int i=1; i<n; i++) {
			xmesh[i] = xmesh[i-1]+dx;
		}
		int nn = data.length;   //nn=N in Botev(2011)
		
		//bin the data uniformly
		int[] hist = histCounts(data, min, max, n);
		
		//scale the bin counts to compute a_k in formula (8.30) via discrete cosine transform.
		double[] initialData = new double[hist.length];
		for (int i=0; i<hist.length; i++) {
			initialData[i] = 1.0*hist[i]/nn;
		}
		double[] a = dct1d(initialData);
		
//		//now compute the optimal bandwidth^2 using the referenced method
//		int[] i2 = new int[n-1];
//		double[] a2 = new double[n-1];
//		for (int i=1; i<n; i++) {
//			i2[i-1] = i*i;
//			a2[i-1] = a[i]*a[i]/4;
//		}
//		//use BracketingNthOrderBrentSolver to solve the equation t=xi*gamma^[7](t) (8.33)
//		//see "root" function defined below
//		double tStar = root(new FixedPoint(nn, i2, a2), nn);
		
		//compute kernel density estimate with parameter tStar via idct
		double tStar = bandwidth * bandwidth / range;
		density = smoothIdct(tStar, a);
		
		//take the rescaling of the data into account
		//since the data was assumed to be in the interval [0,1]
		for (int i=0; i<density.length; i++) {
			if (density[i] > 0) {
				density[i] /= range;
			} else {
				//remove negatives due to round-off error
				density[i] = Math.ulp(1.0);
			}
		}
//		bandwidth = range*Math.sqrt(tStar);
		
//		//cdf estimation
//		//compute optimal bandwidth for cdf estimation
//		double f = fnorm(1, tStar, i2, a2);
//		double tCdf = Math.pow(Math.sqrt(Math.PI)*f*nn, -2.0/3);
//		//obtain density estimate with cdf bandwidth
//		double[] densCdf = smoothIdct(tCdf, a);
//		//now get values of cdf on grid points by summing up values of density
//		cdf = new double[n];
//		cdf[0] = densCdf[0]*0.5*dx/range;
//		for (int i=1; i<n; i++) {
//			cdf[i] = cdf[i-1]+densCdf[i]*dx/range;
//		}
//		//take the rescaling into account
//		bandwidthCdf = range*Math.sqrt(tCdf);
	}

	
	public void kde(double[] data, int n, double min, double max) {
		//round up n to the next power of 2
		n = (int)Math.pow(2, Math.ceil(Math.log(n)/Math.log(2)));
		
		//set up the grid over which the density estimate is computed
		double range = max-min;
		double dx = range/n;
		xmesh = new double[n];
		xmesh[0] = min+dx*0.5;
		for (int i=1; i<n; i++) {
			xmesh[i] = xmesh[i-1]+dx;
		}
		int nn = data.length;//nn=N in Botev(2011)
		
		//bin the data uniformly
		int[] hist = histCounts(data, min, max, n);
		
		//scale the bin counts to compute a_k in formula (8.30) via discrete cosine transf
		double[] initialData = new double[hist.length];
		for (int i=0; i<hist.length; i++) {
			initialData[i] = 1.0*hist[i]/nn;
		}
		double[] a = dct1d(initialData);
		
		//now compute the optimal bandwidth^2 using the referenced method
		int[] i2 = new int[n-1];
		double[] a2 = new double[n-1];
		for (int i=1; i<n; i++) {
			i2[i-1] = i*i;
			a2[i-1] = a[i]*a[i]/4;
		}
		//use BracketingNthOrderBrentSolver to solve the equation t=xi*gamma^[7](t) (8.33)
		//see "root" function defined below
		double tStar = root(new FixedPoint(nn, i2, a2), nn);
		
		//compute kernel density estimate with parameter tStar via idct
		density = smoothIdct(tStar, a);
		
		//take the rescaling of the data into account
		//since the data was assumed to be in the interval [0,1]
		for (int i=0; i<density.length; i++) {
			if (density[i] > 0) {
				density[i] /= range;
			} else {
				//remove negatives due to round-off error
				density[i] = Math.ulp(1.0);
			}
		}
		bandwidth = range*Math.sqrt(tStar);
		
		//cdf estimation
		//compute optimal bandwidth for cdf estimation
		double f = fnorm(1, tStar, i2, a2);
		double tCdf = Math.pow(Math.sqrt(Math.PI)*f*nn, -2.0/3);
		//obtain density estimate with cdf bandwidth
		double[] densCdf = smoothIdct(tCdf, a);
		//now get values of cdf on grid points by summing up values of density
		cdf = new double[n];
		cdf[0] = densCdf[0]*0.5*dx/range;
		for (int i=1; i<n; i++) {
			cdf[i] = cdf[i-1]+densCdf[i]*dx/range;
		}
		//take the rescaling into account
		bandwidthCdf = range*Math.sqrt(tCdf);
	}
	
	private static double[] smoothIdct(double time, double[] a) {
		//smooth the discrete cosine transform of initial data using parameter time
		//in order to apply the inverse discrete cosine transform via equation (8.29)
		double[] at = new double[a.length];
		double s0 = Math.exp(-Math.PI*Math.PI*time/2);
		double s1 = s0;
		double smoother = 1;
		for (int i=0; i<a.length; i++) {
			at[i] = a[i]*smoother;
			smoother *= s1;
			s1 *= s0*s0;
		}
		return idct1d(at);
	}
	
	private static int[] histCounts(double[] data, double min, double max, int numBins) {
		//divides the interval [min,max[ uniformly into numBins bins
		//and counts the number of data entries in each bin
		int[] result = new int[numBins];
		double binSize = (max-min)/numBins;
		for (double d : data) {
			int bin = (int)((d-min)/binSize);
			if ((bin >= 0)&&(bin < numBins)) {
				result[bin] += 1;
			}
		}
		return result;
	}
	
	private static double[] dct1d(double[] data) {
		//computes the discrete cosine transform of the vector data
		//DCT-II with factors adapted according to our kde algorithm (formula 8.30)
		//length of data has to be a power of 2
		int n = data.length;
		//Re-order the elements of the columns of x
		double[] data_reorder = new double[n];
		for (int i=0; i<n/2; i++) {
			data_reorder[i] = data[2*i];
			data_reorder[n/2+i] = data[n-1-2*i];
		}
		Complex[] fft = new FastFourierTransformer(DftNormalization.STANDARD)
				.transform(data_reorder,  TransformType.FORWARD);
		//Multiply FFT by weights
		Complex w0 = new Complex(0,-Math.PI/(2*n)).exp();
		Complex weight = new Complex(2);
		for (int i=1; i<n; i++) {
			weight = weight.multiply(w0);
			fft[i] = fft[i].multiply(weight);
		}
		
		double[] dct = new double[n];
		for (int i=0; i<n; i++) {
			dct[i] = fft[i].getReal();
		}
		return dct;
	}
	
	private static double[] idct1d(double[] data) {
		//computes the inverse discrete cosine transform
		//with factors adapted according to our kde algorithm (formula 8.29)
		//length of data has to be a power of 2
		int n = data.length;
		//Compute weights
		Complex w0 = new Complex(0,Math.PI/(2*n)).exp();
		Complex weight = new Complex(n);
		Complex[] data_weighted = new Complex[n];
		for (int i=0; i<n; i++) {
			data_weighted[i] = weight.multiply(data[i]);
			weight = weight.multiply(w0);
		}
		//Compute x tilde using equation (5.93) in Jain
		Complex[] ifft = new FastFourierTransformer(DftNormalization.STANDARD)
				.transform(data_weighted,  TransformType.INVERSE);
		//Re-order elements of each column according to equations (5.93) and (5.94) in Jain
		double[] idct = new double[n];
		for (int i=0; i<n/2; i++) {
			idct[2*i] = ifft[i].getReal();
			idct[2*i+1] = ifft[n-1-i].getReal();
		}
		return idct;
		//Reference: A. K. Jain, "Fundamentals of Digital Image Processing", pp. 150-153.
	}
	
	private static double fnorm(int s, double time, int[] i2, double[] a2) {
		//this computes an estimate of the squared norm of f^(s) with t_s=time
		//according to formula on p.327
		double f = 0;
		for (int i=0; i<i2.length; i++) {
			f = f + 2*Math.pow(Math.PI, 2*s)*Math.pow(i2[i], s)*a2[i]
					*Math.exp(-i2[i]*Math.PI*Math.PI*time);
		}
		return f;
	}
	
	private static double root(UnivariateFunction f, int nn) {
		//try to find smallest root whenever there is more than one
		if (nn < 50) {
			nn = 50;
		}
		if (nn > 1050) {
			nn = 1050;
		}
		UnivariateSolver solver = new BracketingNthOrderBrentSolver();
		double tol = 1e-12+0.01*(nn-50)/1000;
		boolean flag = false;
		double t = 0;
		do {
			try {
				//try to find a root within the interval [0,tol]
				t = solver.solve(200, f, 0, tol);
				flag = true;
			} catch (NoBracketingException e) {
				tol = tol*2;//double search interval
			} catch (Exception e) {
				tol = 0.1+1e-12;
			}
			if (tol >= 0.1) {
				//if all else fails, find minimum of abs(f)
				UnivariateFunction fAbs = new AbsFunction(f);
				t = new BrentOptimizer(1e-10, 1e-14).optimize(
						new MaxEval(200),
                        new UnivariateObjectiveFunction(fAbs),
                        GoalType.MINIMIZE,
                        new SearchInterval(0,0.1)).getPoint();
				flag = true;
			}
		} while (!flag);
		
		return t;
	}
	
	private static class FixedPoint implements UnivariateFunction {
		//this implements the function t-xi*gamma^[l](t)
		int nn;
		int[] i2;
		double[] a2;
		
		/*input:
		 * nn - size of sample from which density estimate is computed (N in reference)
		 * i2 - array of {1,...,(n-1)^2}
		 * a2 - array of {(a_k/2)^2} where k=1,...,n-1
		 */
		
		FixedPoint(int nn, int[] i2, double[] a2) {
			this.nn = nn;
			this.i2 = i2;
			this.a2 = a2;
		}
		
		public double value(double t) {
			
			int l = 7;
			double f = fnorm(l,t,i2,a2);
			
			for (int s = l-1; s>1; s--) {
				//compute t_s=gamma_s(time) where time=t_{s+1} according to formula (8.32)
				double k0 = 1/Math.sqrt(2*Math.PI);
				for (int k=1; k<2*s; k+=2) {
					k0 *= k;
				}
				double cst = (1+Math.pow(0.5, s+0.5))/3;
				double time = Math.pow(2*cst*k0/nn/f, 2.0/(3+2*s));
				f = fnorm(s,time,i2,a2);
			}
			//we obtained estimate of squared norm of f''
			//now compute t-zeta*t_1=t-xi*gamma^[l](t)
			return t-Math.pow(2*nn*Math.sqrt(Math.PI)*f, -0.4);
		}
	}
	
	private static class AbsFunction implements UnivariateFunction {
		//for a given function f, this implements the function t->abs(f(t))
		UnivariateFunction origFunction;
		
		AbsFunction(UnivariateFunction origFunction) {
			this.origFunction = origFunction;
		}
		
		public double value(double t) {
			return Math.abs(origFunction.value(t));
		}
		
	}
	
}
