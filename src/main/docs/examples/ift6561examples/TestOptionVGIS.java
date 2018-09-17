package ift6561examples;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.probdist.GammaDist;
// import umontreal.ssj.stat.*;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.mcqmctools.*;

/**
 * 
 * @author P. L'Ecuyer
 *
 * European option under VG process, simulated via difference of two gamma variables.
 * Implemements standard MC, IS with exponential twisting, and IS + twisting + conditioning.
 * All with either MC or RQMC.
 */

public class TestOptionVGIS {

    double d = 1;   // One time step.
    double dim = 2;
	double mu = -0.1436;   // Mean BM
	double sigma = 0.12136;  // Volatility of BM
	double theta = mu;   
	double nu;  // Variance parameter of gamma process
	double omega;
	double r;   // Interest rate
	double K;         // Strike
	double s0;        // Initial value of VG
	double T = 1;       // Time horizon
	double discount;    // Discount factor
	double theta0 = 0.0;   // IS twisting parameter
	double muNeg, muPos;  // Means for the two gamma processes 
	double nuNeg, nuPos;  // Var parameters for the two gamma processes 
	double alpha;         // Parameters for the gamma distributions.
	double lambdaNeg, lambdaPos;
	GammaDist distNeg, distPos;  // Gamma distributions under MC 
	GammaDist distNegIS, distPosIS;  // Gamma distributions under twisting
	double factL1;       // Precomputed factor in L for IS1
	double factL2;       // Precomputed factor in L for IS2
    double cis2;         // The constant c for IS with conditional sampling (IS2)

	OptionVG1 vg1 = new OptionVG1 ();
	OptionVGIS1 vgis1 = new OptionVGIS1 ();
	OptionVGIS2 vgis2 = new OptionVGIS2 ();

    public TestOptionVGIS (double mu, double sigma, double nu, double r, double K,
		    double s0, double T, double theta0) {
		this.mu = this.theta = mu;
		this.sigma = sigma;
		this.nu = nu;
		this.r = r;
		this.K = K;
		this.s0 = s0;
		this.T = T;
		this.theta0 = theta0;
		discount = Math.exp(-r * T);
		omega = Math.log(1 - mu * nu - sigma * sigma * 0.5 * nu) / nu;
		muNeg = 0.5 * (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) - theta);
		muPos = 0.5 * (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) + theta);
		nuNeg = muNeg * muNeg * nu;
		nuPos = muPos * muPos * nu;
		alpha = 1.0 / nu;
		lambdaNeg = muNeg / nuNeg;
		lambdaPos = muPos / nuPos;
		distNeg = new GammaDist (alpha, lambdaNeg);
		distPos = new GammaDist (alpha, lambdaPos);
		distNegIS = new GammaDist (alpha, lambdaNeg + theta0);
		distPosIS = new GammaDist (alpha, lambdaPos - theta0);
		factL1 = Math.pow((1.0 - theta0/lambdaPos) * (1.0 + theta0/lambdaNeg), -alpha);
		factL2 = Math.pow((1.0 + theta0/lambdaNeg), -alpha); 
		cis2 = Math.log(K/s0) - r - omega;
	}

	/**
	 * Simulate <tt>m/tt> replications and return the RQMC variance.
	 */
	public static void simulReplicatesRQMCReport (MonteCarloModelDouble model, PointSet p, PointSetRandomization rand,
	        int m, Tally statReps) {
		Chrono timer = new Chrono();
		statReps.init();
		int n = p.getNumPoints();
		// To store the n outputs X for each replication.
		Tally statValue = new Tally();
		PointSetIterator stream = p.iterator();
		for (int rep = 0; rep < m; rep++) {
			rand.randomize(p);
			stream.resetStartStream();
			MonteCarloExperiment.simulateRuns(model, n, stream, statValue);
			statReps.add(statValue.average());   // For the estimator of the mean.
		}
		System.out.println (statReps.report(0.95, 5));
		System.out.printf("Variance per run: %9.5g%n", p.getNumPoints() * statReps.variance());
        System.out.println("Total CPU time:      " + timer.format() + "\n");
	}
	
	public class OptionVG1 implements MonteCarloModelDouble {

		double value;

		/**
		 * Generates the discounted payoff.
		 */
		public void simulate(RandomStream stream) {
			double gneg = distNeg.inverseF(stream.nextDouble());
			double gpos = distPos.inverseF(stream.nextDouble());
			value = discount * (s0 * Math.exp(r + omega + gpos - gneg) - K);
			if (value < 0.0) value = 0.0;
		}

		public double getPerformance() {
			return value;
		}

		public String toString() {
			return "European option on geometric VG process, one time step, MC";
		}
    }
	
	public class OptionVGIS1 extends OptionVG1 {

		public void simulate(RandomStream stream) {
			double gneg = distNegIS.inverseF(stream.nextDouble());
			double gpos = distPosIS.inverseF(stream.nextDouble());
			double L = Math.exp(theta0 * (gneg-gpos)) * factL1;
			value = L * discount * (s0 * Math.exp(r + omega + gpos - gneg) - K); 
			if (value < 0.0) value = 0.0;			
		}
		public String toString() {
			return "European option on geometric VG process, one time step, IS twisting";
		}
	}
	
	public class OptionVGIS2 extends OptionVG1 {

		public void simulate(RandomStream stream) {
			double gneg = distNegIS.inverseF(stream.nextDouble());
			double y = cis2 + gneg;   // We generate gpos conditional on > y
			double umin = distPos.cdf(y); 
			double u = umin + stream.nextDouble() * (1.0-umin);
			double gpos = distPos.inverseF(u);
			double L = factL2 * Math.exp(theta0 * gneg) * 
					   distPos.barF(gneg + cis2);
			value = (L * discount * (s0 * Math.exp(r + omega + gpos - gneg) - K));
			if (value < 0.0) value = 0.0;			
		}
		public String toString() {
			return "European option on geometric VG process, one time step, IS + Cond";
		}
	}
	

	public static void main(String[] args) {
		int dim = 2;
		double T = 1.0;

		double mu = -0.1436;
		double sigma = 0.12136;
		// double theta = mu;
		double r = 0.1;
		double nu = 0.3;
		// double K = 140.0;     // Strike price
		double K = 180.0;     // Strike price
		double s0 = 100.0;
		// double theta0 = 0.0;  // IS twisting parameter   17.32    25.56
        // double theta0 = 17.32;    // for K = 140
        double theta0 = 25.56;    // for K = 180

		TestOptionVGIS test = new TestOptionVGIS (mu, sigma, nu, r, K, s0, T, theta0); 		
		
//		int[] N = { 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288,
//		        1048576, 2097152 }; // 13
//		int[] a = { 149, 275, 857, 1731, 2431, 6915, 12545, 19463, 50673, 96407, 204843, 443165,
//		        768165 };
//		// int[] N = { 16384, 16384, 16384, 16384, 1048576, 1048576, 1048576};
//		// int[] a = { 6229, 6915, 4845, 6063, 401837, 456547, 443165};
//		int mink = 9; // First power of 2 considered.
//		int numSets = 12; // Number of values of n.
//		int numSkipReg = 0;
		int m = 100; // Number of replications.

		// int n = 16 * 1024; // 2^14 for Monte Carlo
		
		RandomStream noise = new MRG32k3a();
		// Chrono timer = new Chrono();
		Tally statValue = new Tally("Stats on value of Asian option for MC");
      

		// SobolSequence pSobol = new SobolSequence(16, 31, dim); // 2^{16} points.


		System.out.println("Pricing an European Option under a geometric VG process.");
		System.out.println ("The exact mean is approximately ");
		System.out.println ("0.10197  for K = 140,  and  1.601e-4  for K = 180");
		System.out.println ("Here we have K = " + K +  ", theta0 = " + theta0);

		// PointSet[] pointSets = new PointSet[numSets];
		// PointSetRandomization rand = new RandomShift(noise);
		// PointSetRandomization randLMS = new LMScrambleShift(new MRG32k3a());
		// RQMCExperiment exper = new RQMCExperiment();
		
		
		// Monte Carlo experiments first.
		int n = 10000000;   // for MC.
		System.out.println("  Ordinary MC:\n");
		MonteCarloExperiment.simulateRunsDefaultReportStudent (test.vg1, n, noise, statValue, 0.95, 4);
		MonteCarloExperiment.simulateRunsDefaultReportStudent (test.vgis1, n, noise, statValue, 0.95, 4);
		MonteCarloExperiment.simulateRunsDefaultReportStudent (test.vgis2, n, noise, statValue, 0.95, 4);

		System.out.println("-----------------------------------------------------\n");

        // Then RQMC experiments.
		System.out.println("  RQMC:\n");
		DigitalNetBase2 pSobol = (new SobolSequence(16, 31, dim-1)).toNetShiftCj();  // 2^{16} points.
		LMScrambleShift randLMS = new LMScrambleShift(new MRG32k3a());
        PointSet pSobolCached = new CachedPointSet(pSobol);

		NestedUniformScrambling randNUS = new NestedUniformScrambling (new MRG32k3a());

		KorobovLattice pKor = new KorobovLattice (65536, 19463, dim);
		BakerTransformedPointSet pKorBaker = new BakerTransformedPointSet(pKor);
		RandomShift randShift = new RandomShift(new MRG32k3a());

		System.out.println ("------------------------------------\n"
				+ "Sobol points with LMS. \n ");
		TestOptionVGIS.simulReplicatesRQMCReport (test.vg1, pSobol, randLMS, m, statValue);
		TestOptionVGIS.simulReplicatesRQMCReport (test.vgis1, pSobol, randLMS, m, statValue);
		TestOptionVGIS.simulReplicatesRQMCReport (test.vgis2, pSobol, randLMS, m, statValue);

		System.out.println ("------------------------------------\n"
				+ "Sobol points with NUS. \n ");
		TestOptionVGIS.simulReplicatesRQMCReport (test.vg1, pSobolCached, randNUS, m, statValue);
		TestOptionVGIS.simulReplicatesRQMCReport (test.vgis1, pSobolCached, randNUS, m, statValue);
		TestOptionVGIS.simulReplicatesRQMCReport (test.vgis2, pSobolCached, randNUS, m, statValue);

		System.out.println ("------------------------------------\n"
				+ "Lattice + bakerS. \n ");
		TestOptionVGIS.simulReplicatesRQMCReport (test.vg1, pKorBaker, randShift, m, statValue);
		TestOptionVGIS.simulReplicatesRQMCReport (test.vgis1, pKorBaker, randShift, m, statValue);
		TestOptionVGIS.simulReplicatesRQMCReport (test.vgis2, pKorBaker, randShift, m, statValue);


		System.out.println("***   THE END   ***\n");

	}
}
