package randvar;

import umontreal.ssj.probdist.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import java.util.Arrays;

public class RandvarExample1 {

	// Generate and print n random variates with generator gen.
	private static void generate(RandomVariateGen gen, int n) {
		double u;
		for (int i = 0; i < n; i++) {
			u = gen.nextDouble();
			System.out.printf("%12.6f%n", u);
		}
		System.out.println("----------------------");
	}

	public static void main(String[] args) {

		// Create four parallel generators (three normal and one gamma).
		RandomVariateGen gen0 = new RandomVariateGen (new MRG31k3p(), new NormalDist());
		RandomVariateGen gen1 = new NormalGen (new MRG31k3p());
		RandomVariateGen gen2 = new NormalGen (new MRG31k3p(), 5.0, 121.4);
		RandomVariateGen gen3 = new GammaGen (new MRG31k3p(), 2.0, 10.0);
		RandomVariateGenInt gen4 = new PoissonGen (new MRG31k3p(), 10.0);
		
		System.out.println ("Some normal, gamma, and Poisson variates \n");
		generate(gen0, 3);  // Generate and prints 3 standard normal variates from gen0
		generate(gen1, 5);  // then 5 more standard normals from gen1
		generate(gen2, 3);  // then 3 normal variates with mean 5 and standard deviation 121.4 
		generate(gen3, 2);  // then 2 gamma variates with parameters (2, 10).
		int[] arrayP = gen4.nextArrayOfInt (10000); // Generate an array of 10000 Poisson variates of mean 10.
		Arrays.sort(arrayP);  
		System.out.println ("50% quantile from Poisson(10) variates: " + arrayP[5000]);		
		System.out.println ("90% quantile from Poisson(10) variates: " + arrayP[9000]);		
		System.out.println ("99% quantile from Poisson(10) variates: " + arrayP[9900]);		
		System.out.println ("99.9% quantile from Poisson(10) variates: " + arrayP[9990]);		
	}
}
