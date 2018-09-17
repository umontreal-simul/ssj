package ift6561examples;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.util.Chrono;

// Extension of Collision class in which the number of items thrown is random 
// with a Poisson distribution with mean mm.
public class CollisionPoisson extends Collision {
	PoissonDist distm;

	public CollisionPoisson(int k, int mm, int maxCounts) {
		super(k, mm, maxCounts);
		distm = new PoissonDist(mm);
	}

	// Here we generate m from the Poisson distribution.
	public int simulate(RandomStream stream) {
		C = 0;
		int loc;
		for (int i = 0; i < k; i++)
			used[i] = false;
		int m = (int) distm.inverseF(stream.nextDouble());
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
		String s = "Collisions in a hashing system\n" + "k = " + k
				+ " locations \n" + "m = " + m
				+ " items on average (Poisson dist.)\n"
				+ "Theorical mean = lambda = " + lambda + "\n";
		return s;
	}

	public static void main(String[] args) {
		int k = 25; // Number of boxes (size of hash table).
		int m = 20;
		int maxCounts = 25; // Values of C >= maxCounts are aggregated.
		int n = 10000 * 1000; // Number of replications.

		Collision col = new Collision(k, m, maxCounts);
		CollisionPoisson colp = new CollisionPoisson(k, m, maxCounts);
		Tally statC = new Tally("Statistics on collisions");
		
		// Run experiment with m fixed.
		System.out.println(col.toString());
		Chrono timer = new Chrono();
		col.simulateRuns(n, new MRG32k3a(), statC);
		System.out.println("Total CPU time:      " + timer.format() + "\n");
		statC.setConfidenceIntervalStudent();
		System.out.println(statC.report(0.95, 3));

		// Run experiment with m Poisson with mean mm.
		timer.init();
		colp.simulateRuns(n, new MRG32k3a(), statC);
		System.out.println("Total CPU time:      " + timer.format() + "\n");
		statC.setConfidenceIntervalStudent();
		System.out.println(statC.report(0.95, 3));

		// Print table with statistics on collision counts and Poisson approximation
		// of P(C=c), for each number c of collisions.
		System.out.println("Counters:\n"
				+ "c  count-fixed  count-poisson  Poisson prob.\n");
		for (int c = 0; c <= maxCounts; c++) {
			System.out.printf("%3d & %8d & %8d & %12.2f  %n", c, col.counts[c],
					colp.counts[c], n * col.poisson.prob(c));
		}
	}
}
