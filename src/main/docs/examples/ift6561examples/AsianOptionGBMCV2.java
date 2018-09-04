package ift6561examples;

import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.stat.Tally;


/**
 * This class represents a <SPAN CLASS="textit">Asian average price call</SPAN> option with European
 * exercise type. It is similar to class <TT>AsianOption</TT>, except that equation is calculated
 * using a <SPAN CLASS="textit">geometric</SPAN> average
 * 
 * <P>
 * </P>
 * <DIV ALIGN="CENTER" CLASS="mathdisplay"> bar(S)<SUB>T</SUB> =
 * (&prod;<SUB>i=1</SUB><SUP>n</SUP><I>S</I>(<I>t</I><SUB>i</SUB>))<SUP>1/n</SUP>. </DIV>
 * <P>
 * </P>
 * 
 * <P>
 * The methods <TT>computeExpectedGeo</TT>, <TT>getPayoffGeo</TT> uses the exact formula from.
 * 
 */
public class AsianOptionGBMCV2 extends AsianOption {

	double[] pathBM;         // Sample path of the underlying BM process.
	double[] muDelta;        // Differences * (r - sigma^2/2)
	double[] sigmaSqrtDelta; // Sqrt(t_{i+1} - t_i)*sigma
	double expectedGeo;    // Expected value of the geometric CV
	double sigma;
	double mu;
	double s0;
	// GeometricBrownianMotion sp;

	/**
	 * Array <TT>obsTimes[0..d+1]</TT> must contain <TT>obsTimes[0]</TT> <SPAN CLASS="MATH">=
	 * 0</SPAN>, plus the <SPAN CLASS="MATH"><I>d</I></SPAN> observation times.
	 * 
	 */
	public AsianOptionGBMCV2(double s0, double r, double sigma, int d, double[] obsTimes,
	        double strike) {

		// Array obsTimes[0..d+1] must contain obsTimes[0]=0, plus the
		// d observation times.
		super(r, d, obsTimes, strike);
		this.s0 = s0;
		this.sigma = sigma;
		mu = r - 0.5 * sigma * sigma;
		muDelta = new double[d + 1];
		sigmaSqrtDelta = new double[d + 1];
		for (int j = 0; j < d; j++) {
			double delta = obsTimes[j + 1] - obsTimes[j];
			muDelta[j] = mu * delta;
			sigmaSqrtDelta[j] = sigma * Math.sqrt(delta);
		}
		computeExpectedGeo();
	}

	/**
	 * Reset the process to <TT>gbm</TT>. Assumes that <TT>obsTimes</TT> have been set.
	 * 
	 */
	public void setProcess(GeometricBrownianMotion gbm) {
		// Reset the process to gbm. Assumes that obsTimes have been set.
		priceProcess = gbm;
		priceProcess.setObservationTimes(obsTimes, d);
	}

	/**
	 * Computes the expected value of geometric CV.
	 * 
	 */
	public double computeExpectedGeo() {
		// Computes the expected value of geometric CV.
		double my = 0;
		double s2y = 0;
		for (int j = 1; j <= d; j++) {
			my += obsTimes[j];
			s2y += (obsTimes[j] - obsTimes[j - 1]) * (d - j + 1) * (d - j + 1);
		}
		my = Math.log(s0) + my * mu / d;
		s2y *= sigma * sigma / (d * d);
		double dd = (-Math.log(strike) + my) / Math.sqrt(s2y);
		return expectedGeo = discount
		        * (Math.exp(my + 0.5 * s2y) * NormalDist.cdf01(dd + Math.sqrt(s2y))
		                - strike * NormalDist.cdf01(dd));
	}

	/**
	 * Returns the expected value of geometric CV.
	 * 
	 */
	public double getExpectedGeo() {
		return expectedGeo;
	}

	/**
	 * Computes and returns discounted payoff for geometric average. Assumes that the path has been
	 * generated.
	 * 
	 */
	public double getPayoffGeo() {
		// Computes and returns discounted payoff for geometric average.
		// Assumes that the path has been generated.
		pathBM = ((GeometricBrownianMotion) priceProcess).getBrownianMotion().getPath();
		double average = 0.0;      // Average over BM sample path.
		for (int j = 1; j <= d; j++)
			average += pathBM[j];
		average /= d;
		average = s0 * Math.exp(average);
		if (average > strike)
			return discount * (average - strike);
		else
			return 0.0;
	}

	/**
	 * Performs <SPAN CLASS="MATH"><I>n</I></SPAN> independent runs using <TT>stream</TT> and
	 * collects statistics in collectors <TT>statValue</TT> and <TT>statValueGeo</TT>.
	 * 
	 */
	public void simulateRuns(int n, RandomStream stream, Tally statValue, Tally statValueGeoCV) {
		// Performs n indep. runs using stream and collects statistics in
		// statValue and statValueGeo.
		statValue.init();
		statValueGeoCV.init();
		for (int i = 0; i < n; i++) {
			path = priceProcess.generatePath(stream);
			statValue.add(getPerformance());
			statValueGeoCV.add(getPayoffGeo() - expectedGeo);
			stream.resetNextSubstream();
		}
	}

	/**
	 * Makes <SPAN CLASS="MATH"><I>m</I></SPAN> independent randomizations of the point set
	 * <SPAN CLASS="MATH"><I>p</I></SPAN> using stream <TT>noise</TT>. For each of them, performs
	 * one simulation run for each point of <SPAN CLASS="MATH"><I>p</I></SPAN>, and adds the
	 * averages over these points to the collectors <TT>statRQMC</TT> and <TT>statRQMCCV</TT>.
	 * 
	 */
	public void simulateRQMCCV(int m, PointSet p, RandomStream noise, Tally statRQMC,
	        Tally statRQMCCV) {
		// Makes m independent randomizations of the point set p using stream
		// noise. For each of them, performs one simulation run for each point
		// of p, and adds the averages over these points to the collectors
		// statRQMC and statRQMCCV.

		Tally statValue = new Tally("stat on value");
		Tally statValueGeoCV = new Tally("stat on value of CV");
		statRQMC.init();
		statRQMCCV.init();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			p.randomize(noise);
			// p.addRandomShift (0, d, noise);
			stream.resetStartStream();
			simulateRuns(p.getNumPoints(), stream, statValue, statValueGeoCV);
			statRQMC.add(statValue.average());
			statRQMCCV.add(statValueGeoCV.average());
		}
	}

}
