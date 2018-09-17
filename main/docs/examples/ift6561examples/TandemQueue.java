package ift6561examples;
/* Simulates a tandem queue as in the class notes of ift6561.
 The constructor takes both a time limit T and a max number of customers N_c,
 and the simulation stops whenever one of these limits is reached.
 Stations 0 (which represents the arrival process) and station 1 must have 
 infinite capacity. The vectors mu and capacity must have dimension m+1.
 The matrix D is stored in circular buffers (arrays) of just the right sizes
 to store what is needed. Two indexes give pointers to where are the values of 
 D_{j,i} for the current and previous customers, in the buffer, for each station. 
 */

import umontreal.ssj.probdist.ExponentialDist;
import umontreal.ssj.charts.HistogramChart;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.TallyStore;
import umontreal.ssj.rng.MRG32k3a;

public class TandemQueue {

    // System parameters.
    double lambda;    // Rate of interarrival times distribution.
    int m;    	      // Number of service stations.
    double[] mu;      // Service rate at each station; mu[j] for station j.
    int[] capacity;   // Capacity at each station; capacity[j] for station j.
                      // capacity[j] = 0 means an infinite capacity.

    // Working variables.
    double[][] D;	// Circular buffers to store the arrival times D_{0,i}
                        // and the departure times D_{j,i} from station j.
    int[] sizeD;        // sizeD[j] is the size of D[j] for each j;
                        // = capacity[j] if finite, = 1 otherwise.
    int[] id;		// Index of current customer i in circular array D[j].
    int[] idprev;	// Index of previous customer i-1 in D[j].
    double[] W;		// Cumulated total waiting time at each station.
    double[] B;	        // Cumulated total blockage time at each station.
    RandomStream streamA;  // Stream for interarrival times.
    RandomStream streamS;  // Stream for service times.

    public TandemQueue(double lambda, int m, double[] mu, int[] capacity) {
        this.lambda = lambda;
        this.m = m;
        this.mu = mu;
        this.capacity = capacity;
        capacity[0] = capacity[1] = 0; // Force infinite capacities there. 
        sizeD = new int[m + 1];
        W = new double[m + 1];
        B = new double[m + 1];
        // We will store the last sizeD[j] departure times for each station j.
        D = new double[m + 1][];
        for (int j = 0; j <= m; j++) {
            if (capacity[j] == 0) {
                sizeD[j] = 1;
            } else {
                sizeD[j] = capacity[j];
            }
            D[j] = new double[sizeD[j]];
        }
        id = new int[m + 1];      // Index of current i in D[j] is id[j].
        idprev = new int[m + 1];  // Index of previous i in D[j] is idprev[j].
        streamA = new MRG32k3a();
        streamS = new MRG32k3a();
    }
 
    // Generates the trajectory of one customer through the system.
    public void oneCustomerPath() {
        double sji;	// Service time.
        double wji;	// Waiting time.
        double bji;     // Blockage time.
        for (int j = 1; j <= m; j++) {
            // Compute wji first, because D[j][idprev[j]] will change if size[j]=1.
            wji = Math.max(0, D[j][idprev[j]] - D[j - 1][id[j - 1]]);
            sji = ExponentialDist.inverseF(mu[j], streamS.nextDouble());
            
            // We first compute the departure time of Cust. i (whose index at 
            // station j is id[j]) without taking the blocking time into account.
            D[j][id[j]] = sji + Math.max(D[j - 1][id[j - 1]], D[j][idprev[j]]);
            if (j < m) {
                // If the next station (j+1) has finite capacity, we look for the
            	// departure time of cust. i-c_{j+1} at that station.  
            	// The index of customer i-c_{j+1} in the circular table for station j+1 
            	// is id[j+1], the same as cust. i, so the departure time we look for is
            	// D[j+1][id[j+1]].  If it is larger than D[j][id[j]],
            	// we update  D[j][id[j]] to this value.
                // D[j + 1][id[j + 1]] will be updated at the next iteration.
                if (capacity[j + 1] > 0) {
                    D[j][id[j]] = Math.max(D[j][id[j]], D[j + 1][id[j + 1]]);                  
                }
            }
            bji = D[j][id[j]] - D[j - 1][id[j - 1]] - wji - sji;
            W[j] += wji;
            B[j] += bji;
        }
    }

    //Simulates the system once.
    public void simulateOneRun(int maxNc, double maxT) {
        // Initialize counters.
        for (int j = 0; j <= m; j++) {
            id[j] = idprev[j] = 0;
            for (int i = 0; i < sizeD[j]; i++) {
                D[j][i] = 0.0;
            }
            W[j] = B[j] = 0.0;
        }
        double Ti = ExponentialDist.inverseF(lambda, streamA.nextDouble());
        for (int i = 1; (i < maxNc) & (Ti < maxT); i++) {
            oneCustomerPath();
            D[0][0] = Ti += ExponentialDist.inverseF(lambda, streamA.nextDouble());
            // Update the index for departure times.
            for (int j = 2; j <= m; j++) {
               idprev[j] = id[j];
               id[j]++;
               if (id[j] >= sizeD[j]) {
                  id[j] = 0;
               }
            }
        }
    }

    // Simulate the system n times.
    public void simulateRuns(int n, int maxNc, double maxT) {
        TallyStore wait[] = new TallyStore[m + 1];
        TallyStore block[] = new TallyStore[m + 1];
        for (int j = 1; j <= m; j++) {
            wait[j] = new TallyStore("Total waiting time at station " + j);
            block[j] = new TallyStore("Total blocking time at station " + j);
        }
        for (int i = 0; i < n; i++) {
            simulateOneRun(maxNc, maxT);
            for (int j = 1; j <= m; j++) {
                wait[j].add(W[j]);
                block[j].add(B[j]);
            }
        }
        displayResults(wait, m);
        displayResults(block, m - 1);
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
        int capacity[] = {0, 0, 4, 8};
        int maxNc = 1000000;
        double maxT = 1000.0;    // Fixed time horizon.
        TandemQueue tandem = new TandemQueue(lambda, m, mu, capacity);
        tandem.simulateRuns(1000, maxNc, maxT);
    }
}