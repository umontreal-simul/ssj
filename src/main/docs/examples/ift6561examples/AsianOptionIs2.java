package ift6561examples;
// Exercice 2 (e)
// package travail_final;
import umontreal.ssj.probdist.ContinuousDistribution;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.TruncatedDist;
import umontreal.ssj.randvar.GammaAcceptanceRejectionGen;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class AsianOptionIs2 {

	double strike; // Strike price.
	int d; // Number of observation times.
	double r; // taux d inetret
	double sigma;// volatilit� du MB
	double theta;// moyenne du MB
	double s0;
	double nu;// Variance du processus Gamma
	double T;// horizon de temps
	double theta_para;
	double discount;

	public AsianOptionIs2(double nu, double theta, double r, double sigma,
			double strike, double s0, int d, double T, double theta_para) {
		this.strike = strike;
		this.d = d;
		this.sigma = sigma;
		this.theta = theta;
		this.T = T;
		this.s0 = s0;
		discount = Math.exp(-r * T);
		this.nu = nu;
		this.r = r;
		this.theta_para = theta_para;
	}

	public void simulateIs(RandomStream stream, Tally statIS, double[] z, int d) {
		double mup = (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) + theta) / 2;
		double mun = (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) - theta) / 2;
		double mupg = (1 / nu);
		double volpg = (1 / (mup * nu));
		double mung = (1 / nu);
		double volng = (1 / (mun * nu));

		// On commence par estime en utilisant monte Carlo Standard
		double w = Math.log(1 - theta * nu - sigma * sigma * 0.5 * nu) / nu;
		// System.out.println(w);
		double s = Math.log(strike / s0);
		GammaAcceptanceRejectionGen gnegIS = new GammaAcceptanceRejectionGen(
				stream, new GammaDist(mung, volng + theta_para));
		// System.out.println(stream.nextDouble());
		double y1IS = gnegIS.nextDouble();// On prend GIS+(1)
		ContinuousDistribution dist = new GammaDist(mupg, volpg);
		TruncatedDist dist2 = new TruncatedDist(dist, y1IS + s - r - w, 100000);
		double y2IS = dist2.inverseF(stream.nextDouble());
		// Rapport de vraisemblance
		double L = Math.exp(theta_para * (y1IS))
				* Math.pow((volng / (volng + theta_para)), mung)
				* dist.barF(y1IS + s - r - w);

		statIS.add(L * discount * (s0 * Math.exp(r + w + y2IS - y1IS) - strike));

	}

	public void simulateISNRUNS(int N, RandomStream stream, Tally statIS,
			double[] z, int d) {
		statIS.init();

		for (int l = 1; l < N; l++) {
			simulateIs(stream, statIS, z, d);
			stream.resetNextSubstream();
		}

	}

	// simulation par MC standard
	public void simulateMcNRUNS(int N, RandomStream stream, Tally statMC,
			double[] z, int d) {

		statMC.init();
		double mup = (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) + theta) / 2;
		double mun = (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) - theta) / 2;

		GammaAcceptanceRejectionGen gpos = new GammaAcceptanceRejectionGen(
				stream, new GammaDist(1 / nu, 1 / (mup * nu)));
		GammaAcceptanceRejectionGen gneg = new GammaAcceptanceRejectionGen(
				stream, new GammaDist(1 / nu, 1 / (mun * nu)));

		for (int l = 1; l <= N; l++) {
			// La m�thode monte carlo standard
			double w = Math.log(1 - theta * nu - sigma * sigma * 0.5 * nu) / nu;
			// System.out.println(w);
			double s = Math.log(strike / s0);
			double y1 = gneg.nextDouble();// On prend G+(1)
			double y2 = gpos.nextDouble();

			if (s <= r + w + y2 - y1)
				statMC.add(discount * (s0 * Math.exp(r + w + y2 - y1) - strike));
			else
				statMC.add(0.0);

		}
	}

	public static void main(String[] args) {
		int N = 100000000;
		int d = 1;
		double[] z = new double[d + 1];

		for (int t = 1; t <= d; t++)
			z[t] = (double) t / (double) d;

		double r = 0.1;
		double sigma = 0.12136;
		double theta = -0.1436;
		double s0 = 100;
		double nu = 0.3;
		double T = 1;

		Tally statMC = new Tally(
				"Stats of the the estimation with standard MC ");
		Tally statIS = new Tally("Stats of the the estimation with IS ");
		// k=180,theta= 25.5652
		// k=140,THETA=17.3194

		AsianOptionIs2 Ais2 = new AsianOptionIs2(nu, theta, r, sigma, 180,
				s0, d, T, 25.5652);
		Ais2.simulateISNRUNS(N, new MRG32k3a(), statIS, z, d);
		System.out.println("Pour une valeur du strike �gale � " + 180);
		statIS.setConfidenceIntervalStudent();
		System.out.println(statIS.report(0.95, 3));
		Ais2.simulateMcNRUNS(N, new MRG32k3a(), statMC, z, d);
		statMC.setConfidenceIntervalStudent();
		System.out.println(statMC.report(0.95, 3));
		System.out
				.println(" Le rapport de r�duction de la variance avec la m�thode IS est "
						+ statMC.variance() / statIS.variance());

	}
}
