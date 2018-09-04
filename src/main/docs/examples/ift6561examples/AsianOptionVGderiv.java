package ift6561examples;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.rng.*;


public class AsianOptionVGderiv {

	StochasticProcess priceProcess; // Underlying process for the price.
	int d; // Number of observation times.
	double[] obsTimes; // Observation times.
	double[] path; // Sample path of the process.
	double strike; // Strike price.
	double discount; // Discount factor exp(-r * obsTimes[t]).

	// Array obsTimes[0..d+1] must contain obsTimes[0]=0.0,
	// plus the d observation times.

	public AsianOptionVGderiv (double r, int d, double[] obsTimes, double strike) {
		this.d = d;
		this.obsTimes = new double[d + 1];
		for (int j = 0; j <= d; j++)
			this.obsTimes[j] = obsTimes[j];
		this.strike = strike;
		discount = Math.exp(-r * obsTimes[d]);
	}

	
	public AsianOptionVGderiv (StochasticProcess sp, double r, int d,
			double[] obsTimes, double strike) {
		this(r, d, obsTimes, strike);
		setProcess(sp);
	}


	public void setProcess(StochasticProcess sp) {
		// Reset the process to sp. Assumes that obsTimes have been set.
		priceProcess = sp;
		sp.setObservationTimes(obsTimes, d);
	}

	/**
	 * Computes and returns discounted payoff. Assumes path has been generated.
	 * 
	 */
	public double getValue() {
		// Computes and returns discounted payoff. Assumes path has been
		// generated.
		double average = 0.0; // Average over sample path.
		for (int j = 1; j <= d; j++)
			average += path[j];
		average /= d;
		if (average > strike)
			return discount * (average - strike);
		else
			return 0.0;
	}



	public void simulate(RandomStream stream) {
		path = priceProcess.generatePath(stream);
	}

	// Performs n indep. runs using stream and collects statistics.
	public void simulateRuns(int n, RandomStream stream, Tally statValue) {
		statValue.init();
		for (int i = 0; i < n; i++) {
			simulate(stream);
			statValue.add(getValue());
			stream.resetNextSubstream();
		}
	}

	// Performs n indep. runs using stream and collects statistics.
	// The second collector collects only the positive payoffs.
	public void simulateRuns(int n, RandomStream stream, Tally statValue,
			Tally statValuePos) {
		statValue.init();
		statValuePos.init();
		double x;
		for (int i = 0; i < n; i++) {
			simulate(stream);
			x = getValue();
			statValue.add(x);
			if (x > 0.0000000001)
				statValuePos.add(x);
			stream.resetNextSubstream();
		}
	}

	public String toString() {
		return "Asian option model with " + d + " observation times";
	}

	public static void main(String[] args) {
		int d = 16;
		double[] obsTimes = new double[d + 1];
		obsTimes[0] = 0.0;
		for (int j = 1; j <= d; j++)
			obsTimes[j] = (double) j / (double) d;		
		int n = 100000;
		AsianOptionVGderiv asian = new AsianOptionVGderiv (0.1, d, obsTimes, 101.0);
		RandomStream stream = new MRG32k3a();
		double mu = -0.1436; 			
		double sigma = 0.12136; 		
		double theta = mu;
		double r = 0.1; 				
		double s0 = 100.0; 				
		double v = 0.3;	
		
		
		// BGSS
		
		//StochasticProcess process = new GeometricVarianceGammaProcess(s0, r,
		//		new VarianceGammaProcess(s0, theta, sigma, v, stream));
		
		// BGBS
		
		StochasticProcess process = new GeometricVarianceGammaProcess(s0, r,
				new VarianceGammaProcess(s0, new BrownianMotionBridge(0,mu,sigma,stream),
						new GammaProcessSymmetricalBridge(s0,1,v,stream)));
		// DGBS
		  
		/*StochasticProcess process = new GeometricVarianceGammaProcess(s0, r,
				new VarianceGammaProcessDiffPCABridge(s0, theta, sigma, v, stream));
		*/
		
		
		asian.setProcess(process);
		Tally statValue = new Tally("Stats on value of Asian option");
		asian.simulateRuns(n, stream, statValue);
		statValue.setConfidenceIntervalStudent();
		System.out.println(statValue.report(0.95, 4));
		System.out.printf("Var. per run: %9.4g%n", statValue.variance());
		
		//Derivative for s0
				stream.resetStartStream(); //Use common numbers
				double delta_s0 = 0.000001; 
				Tally statValueS0 = new Tally("Stats on value of Asian option: s0+epsi_s0 ");
				s0 = s0 + delta_s0;
				process = new GeometricVarianceGammaProcess(s0, r,
						new VarianceGammaProcess(s0, new BrownianMotionBridge(0,mu,sigma,stream),
								new GammaProcessSymmetricalBridge(s0,1,v,stream)));
			//	process = new GeometricVarianceGammaProcess(s0, r,
			//			new VarianceGammaProcess(s0, theta, sigma, v, stream));
				asian.setProcess(process);
				asian.simulateRuns(n, stream, statValueS0);
				System.out.println("Derivative for S0");
				System.out.println((statValueS0.average() - statValue.average())/delta_s0);
				
		//Derivative for v
				stream.resetStartStream();
				double delta_v  = 0.000001;
				Tally statValueV = new Tally("Stats on value of Asian option: v+delta_v ");
				s0 = s0 - delta_s0; // Go back to initial s0
				v = v + delta_v;
				 process = new GeometricVarianceGammaProcess(s0, r,
						new VarianceGammaProcess(s0, new BrownianMotionBridge(0,mu,sigma,stream),
								new GammaProcessSymmetricalBridge(s0,1,v,stream)));
				//process = new GeometricVarianceGammaProcess(s0, r,
				//		new VarianceGammaProcess(s0, theta, sigma, v, stream));
				asian.setProcess(process);
				asian.simulateRuns(n, stream, statValueV);
				System.out.println("Derivative for v");
				System.out.println((statValueV.average() - statValue.average())/delta_v);
				
		
	}
}
