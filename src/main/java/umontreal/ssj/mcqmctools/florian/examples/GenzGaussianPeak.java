package umontreal.ssj.mcqmctools.florian.examples;

import java.util.Arrays;

import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.rng.RandomStream;
/**
 * This class implements the Gaussian peak function of the Genz function package,
 * see \link https://www.sfu.ca/~ssurjano/gaussian.html \endlink. This function is given by
 * 
 * \f[ f(x_1,x_2,\dots,x_d) = \exp\left( - \sum_{j=1}^d a_j^2 (x_j-u_j)^2 \right) \f],
 * 
 * where @f$x_j, u_j@f\in[0,1]$.
 * @author florian
 *
 */
public class GenzGaussianPeak implements MonteCarloModelDouble{

	double dim;
	double[] a;
	double[] u;
	double performance;
	
	/**
	 * Constructor passing the shape parameters \a a, \a u, and @f$d@f$.
	 * @param a
	 * @param u
	 * @param dim the value for @f$d@f$.
	 */
	public GenzGaussianPeak(double[] a, double[] u,int dim) {
		this.a = new double [dim];
		this.u = new double [dim];
		for(int j =0; j < dim; ++j) {
			this.a[j] = a[j]; 
			this.u[j] = u[j]; 
		}
		this.dim = dim;
	}
	
	/**
	 * Constructor taking @f$d@f$ as the minimum of the length of a and of u.
	 * @param a
	 * @param u
	 */
	public GenzGaussianPeak(double[] a,double[] u) {
		this(a,u, Math.min(a.length,u.length));
	}
	
	/**
	 * Constructor taking @f$a_j=a@f$ and @f$u_j=u@f$ for all @f$1\leq j\leq d@f$.
	 * @param a
	 * @param u
	 * @param dim the value for @f$d@f$.
	 */
	public GenzGaussianPeak(double a, double u,int dim) {
		this.a = new double[dim];
		Arrays.fill(this.a,a);
		this.u = new double[dim];
		Arrays.fill(this.u,u);
		this.dim = dim;
	}
	
	public void simulate(RandomStream stream) {
		 double exponent = 0.0;
		for(int j =0; j < dim; ++j)
			exponent += a[j]*a[j]*Math.pow(stream.nextDouble() - u[j],2);
		performance = Math.exp(-exponent);
	}

	public double getPerformance() {
		return performance;
	}
	
	public String toString() {
		String str = "GenzGaussianPeak [a = {" + a[0];
		for(int j = 1; j < dim; ++j)
			str += ", " +a[j];
		str+="}, ";
		str+= "u = {" + u[0];
				for(int j = 1; j < dim; ++j)
					str += ", " +u[j];
		str+="}]\n";
		return str;
	}

}
