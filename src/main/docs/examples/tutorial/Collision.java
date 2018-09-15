package tutorial;

import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.Chrono;

public class Collision {
   int k;            // Number of locations.
   int m;            // Number of items.
   double lambda;    // Theoretical expectation of C (asymptotic).
   boolean[] used;   // Locations already used.
   int maxCounts; // Values of C >= maxCounts are aggregated.
   int[] counts; // Counts the number of occurrences of each value of C.
   PoissonDist poisson; // Will be a Poisson distribution with mean lambda.
	
   public Collision (int k, int m, int maxCounts) {
      this.k = k;
      this.m = m;
      lambda = (double) m * m / (2.0 * k);
      used = new boolean[k];
	  this.maxCounts = maxCounts;
	  counts = new int[maxCounts + 1];
	  poisson = new PoissonDist(lambda);
   }

   // Generates and returns the number of collisions.
   public int simulate (RandomStream stream) {
      int C = 0;
      for (int i = 0; i < k; i++) used[i] = false;
      for (int j = 0; j < m; j++) {
         int loc = stream.nextInt (0, k-1);
         if (used[loc]) C++;
         else used[loc] = true;
      }
      return C;
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

   public static void main (String[] args) {
	    int k = 10000;  int m = 500;
		int maxCounts = 30;
		int n = 10000000;
		Collision col = new Collision(k, m, maxCounts);
		Tally statC = new Tally("Statistics on collisions");
		// System.out.println(col.toString());
		Chrono timer = new Chrono();
		col.simulateRuns(n, new MRG32k3a(), statC);
		System.out.println("Total CPU time:      " + timer.format() + "\n");
		statC.setConfidenceIntervalStudent();
		System.out.println(statC.report(0.95, 3));
		System.out.println ("Theoretical mean: lambda = " + col.lambda + "\n");

		System.out.println("Counters:\n"
				+ "  c       count    Poisson expect.\n");
		for (int c = 0; c <= col.maxCounts; c++) {
			System.out.printf(" %2d  %10d   %12.1f%n", c, col.counts[c],
					n * col.poisson.prob(c));
		}
	}
}
