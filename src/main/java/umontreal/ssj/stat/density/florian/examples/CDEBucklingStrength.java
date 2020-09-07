package umontreal.ssj.stat.density.florian.examples;

import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.density.ConditionalDensityEstimator;

/**
 * This class implements a conditional density estimator for the buckling
 * strength @f$X@f$ of a steel plate as given by Schields and Zhang '16. The
 * model is given by
 * 
 * \f[ X=\left( \frac{2.1}{\lambda}-\frac{0.9}{\lambda^2} \right) \left(
 * 1-\frac{0.75\delta_0}{\lambda}\right) \left(1-\frac{2\eta t}{b} \right), \f]
 * 
 * 
 * where \lambda = (b/t) \sqrt{\sigma_0/eta} and where the
 * parameters @f$b,t,\sigma_0,E,\delta,\eta@f$ are random variables. More
 * precisely, $t$ and $b_0$ follow a lognormal distribution and all the others
 * are normally distributed.
 * 
 * The realizations of the random variables at which this CDE can be constructed can 
 * be obtained with the class #BucklingStrengthVars.
 * 
 * This estimator consists of two individual estimators @f$\hat{f}_{\delta_0}@f$
 * and @f$\hat{f}_{\eta}@f$, where the index indicates which variable is hidden.
 * For the final estimator we take a convex combination of the above estimators, i.e.,
 * @f$\hat{f} = p \hat{f}_{\delta_0}+(1-p)\hat{f}_{\eta}@f$.
 * 
 * 
 * 
 * @author florian
 *
 */

public class CDEBucklingStrength extends ConditionalDensityEstimator {

	double mud, mue;
	double sigmad, sigmae;
	double p;

	/**
	 * Constructor
	 * @param mud mean of @f$\delta_0@f$.
	 * @param sigmad standard deviation of @f$\delta_0@f$
	 * @param mue mean of @f$\eta@f$.
	 * @param sigmae standard deviation of @f$\eta@f$.
	 * @param p factor for convex combination.
	 */
	public CDEBucklingStrength(double mud, double sigmad, double mue, double sigmae, double p) {
		this.mue = mue;
		this.sigmae = sigmae;
		this.mud = mud;
		this.sigmad = sigmad;
		this.p = p;
	}

	@Override
	public double evalEstimator(double x, double[] data) {

		double fac1 = 2.1 / data[1] - 0.9 / (data[1] * data[1]);
		double fac2 = 1.0 - data[0] * data[3];
		double fac3 = 1.0 - 0.75 * data[2] / data[1];
		double fac4 = data[1] / (0.75);
		double fac5 = 1.0 / (data[0]);

		if (fac1 < 0 || fac2 < 0 || fac3 < 0 || fac4 < 0 || fac5 < 0) {
//			System.out.println(fac1 + "\t" + fac2 + "\t" + fac3 + "\t" + fac4 + "\t" + fac5);
			return 0.0;
		}

		else
			return (p * NormalDist.density(mud, sigmad, (1.0 - x / (fac1 * fac2)) * fac4) * fac4 / (fac1 * fac2)
					+ (1.0 - p) * NormalDist.density(mue, sigmae, (1.0 - x / (fac1 * fac3)) * fac5) * fac5
							/ (fac1 * fac3));
	}

	public String toString() {
		return "CDEBucklingStrength";
	}

	public double getP() {
		return p;
	}

	public void setP(double p) {
		this.p = p;
	}

}
