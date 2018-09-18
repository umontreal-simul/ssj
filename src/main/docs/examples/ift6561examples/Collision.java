package ift6561examples;

import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.util.Chrono;

// We throw m balls at random in k locations and count 
// the number of C collisions.  We repeat n times and count the
// frequency of each value of C.

public class Collision {
	int k; // Number of locations.
	int m; // Number of items (balls) that are thrown.
	int C; // Number of collisions.
	double lambda; // Theoretical expectation of C (asymptotic approx.).
	boolean[] used; // Locations already used.
	int maxCounts; // Values of C >= maxCounts are aggregated.
	int[] counts; // Counts the number of occurrences of each value of C.
	PoissonDist poisson; // Will be a Poisson distribution with mean lambda.

	public Collision(int k, int m, int maxCounts) {
		this.k = k;
		this.m = m;
		lambda = (double) m * m / (2.0 * k);
		used = new boolean[k];
		this.maxCounts = maxCounts;
		counts = new int[maxCounts + 1];
		poisson = new PoissonDist(lambda);
	}

	// Generates the number of collisions, for fixed m.
	public int simulate(RandomStream stream) {
		C = 0;
		int loc;
		for (int i = 0; i < k; i++)
			used[i] = false;
		for (int j = 0; j < m; j++) {
			loc = stream.nextInt(0, k - 1);
			if (used[loc] == true)
				C++;
			else
				used[loc] = true;
		}
		return C;
	}

	public String toString() {
		String s = "Collisions in a hashing system\n" + "k = " + k + " locations \n" + "m = " + m
		        + " items \n" + "Theorical mean = lambda = " + lambda + "\n";
		return s;
	}

	// Performs n indep. runs using stream and collects statistics in statC.
	public void simulateRuns(int n, RandomStream stream, Tally statC) {
		statC.init();
		int C;
		for (int c = 0; c < maxCounts; c++)
			counts[c] = 0;
		for (int i = 0; i < n; i++) {
			C = simulate(stream);
			statC.add(C);
			if (C > maxCounts)
				C = maxCounts;
			counts[C]++;
		}
	}

	public static void main(String[] args) {
		int k = 100;
		int m = 40;
		int maxCounts = 30;
		int n = 10000000;
		Collision col = new Collision(k, m, maxCounts);
		Tally statC = new Tally("Statistics on collisions");
		System.out.println(col.toString());
		Chrono timer = new Chrono();
		col.simulateRuns(n, new MRG32k3a(), statC);
		System.out.println("Total CPU time:      " + timer.format() + "\n");
		statC.setConfidenceIntervalStudent();
		System.out.println(statC.report(0.95, 3));

		System.out.println("Counters:\n" + "c  count  fraction  Poisson prob.\n");
		for (int c = 0; c <= col.maxCounts; c++) {
			System.out.printf(c + "  " + col.counts[c] + "  %10.6g%n", n * col.poisson.prob(c));
		}
	}
}
