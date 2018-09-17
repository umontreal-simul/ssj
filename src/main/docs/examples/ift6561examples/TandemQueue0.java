package ift6561examples;
/* Simulates a tandem queue as in the notes of ift6561, with no blocking.
 The constructor takes both a time limit T and a max number of customers N_c,
 and the simulation stops whenever one of these limits is reached.  
*/

import umontreal.ssj.probdist.ExponentialDist;
import umontreal.ssj.charts.HistogramChart;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.rng.MRG32k3a;

public class TandemQueue0 {

    // System parameters.
    double lambda;    // Rate of interarrival times distribution.
    int m;    	      // Number of service stations.
    double[] mu;      // Service rate at each station; mu[j] for station j.

    // Working variables.
    double[] D;	      // Stores the departure times D_{j,i} from station j.
    double[] W;		// Cumulated total waiting time at each station.
    double[] B;	        // Cumulated total blockage time at each station.
    RandomStream streamA;  // Stream for interarrival times.
    RandomStream streamS;  // Stream for service times.

    public TandemQueue0(double lambda, int m, double[] mu) {
        this.lambda = lambda;
        this.m = m;
        this.mu = mu;
        W = new double[m + 1];
        B = new double[m + 1];
        D = new double[m + 1];
        streamA = new MRG32k3a();
        streamS = new MRG32k3a();
    }
 
    // Generates the trajectory of one customer through the system.
    public void oneCustomerPath() {
        double sji;	// Service time.
        double wji;	// Waiting time.
        for (int j = 1; j <= m; j++) {
            // Compute wji first, because D[j] will change.
            wji = Math.max(0, D[j] - D[j - 1]);
            W[j] += wji;
            sji = ExponentialDist.inverseF(mu[j], streamS.nextDouble());
            D[j] = sji + Math.max(D[j - 1], D[j]);
        }
    }

    //Simulates the system once.
    public void simulateOneRun(int maxNc, double maxT) {
        // Initialize counters.
        for (int j = 0; j <= m; j++) {
            D[j] = W[j] = B[j] = 0.0;
        }
        double Ti = ExponentialDist.inverseF(lambda, streamA.nextDouble());
        for (int i = 1; (i < maxNc) & (Ti < maxT); i++) {
            oneCustomerPath();
            D[0] = Ti += ExponentialDist.inverseF(lambda, streamA.nextDouble());
        }
    }

    // Simulate the system n times.
    public void simulateRuns(int n, int maxNc, double maxT) {
        TallyStore wait[] = new TallyStore[m + 1];
        for (int j = 1; j <= m; j++) {
            wait[j] = new TallyStore("Total waiting time at station " + j);
        }
        for (int i = 0; i < n; i++) {
            simulateOneRun(maxNc, maxT);
            for (int j = 1; j <= m; j++) {
                wait[j].add(W[j]);
            }
        }
        displayResults(wait, m);
    }

    // Display the results. 
    public void displayResults(TallyStore[] tally, int range) {
        for (int j = 1; j <= range; j++) {
            tally[j].setConfidenceIntervalStudent();
            System.out.println(tally[j].report(0.95, 3));
            HistogramChart chart;
            double[] data = tally[j].getArray();
            chart = new HistogramChart(tally[j].getName(), null, null, data);
            chart.view(800, 500);
            chart.toLatexFile(tally[j].getName() + ".tex", 12, 8);

        }
    }

    public static void main(String[] args) {
        int m = 3;
        double lambda = 1.0;
        double mu[] = {0.0, 1.5, 1.2, 1.2};
        int maxNc = 1000000;
        double maxT = 1000.0;    // Fixed time horizon.
        TandemQueue0 tandem = new TandemQueue0(lambda, m, mu);
        tandem.simulateRuns(10000, maxNc, maxT);
    }
}
