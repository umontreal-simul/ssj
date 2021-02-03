package umontreal.ssj.markovchainrqmc.florian.chemicalReactionNetworks;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.markovchainrqmc.florian.chemicalReactionNetworks.ChemicalReactionNetwork;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.util.Chrono;

import umontreal.ssj.util.sort.MultiDim01;
import umontreal.ssj.util.sort.florian.MultiDim;
/**
 * Implements the Schloegl system, see, e.g., Beentjes and Baker '18. Since the total number of molecules does not change through any
 * reaction, the state space is, in fact, two dimensional.
 * @author florian
 *
 */
public class SchloeglSystem extends ChemicalReactionNetwork implements MultiDim{

	/**
	 * total number of molecules
	 */
	double N0; // Total number of molecules
	/**
	 * Means of the states over 15 steps for the example in our paper
	 */
	double[][] means = {
			{0.0,0.0},
			{250.388, 99819.6},
			{250.876, 99638.8},
			{251.564, 99457.4},
			{252.604, 99275.1},
			{253.987, 99091.4},
			{255.527, 98906.4},
			{256.891, 98720.2},
			{257.754, 98533.2},
			{257.838, 98346.3},
			{257.059, 98160.2},
			{255.55, 97975.6},
			{253.3, 97793.2},
			{250.431, 97613.6},
			{247.077, 97437.2}
	}
	; 
	/**
	 * Standard deviations of the states over 15 steps for the example in our paper
	 */
	double [][] stdDevs = {
			{1.0,1.0},
			{24.5894, 17.8271},
			{38.9162, 46.7408},
			{53.3392, 89.9386},
			{68.5439, 147.717},
			{84.3616, 220.439},
			{100.014, 307.709},
			{114.611, 408.377},
			{127.464, 520.577},
			{138.278, 642.256},
			{146.982, 771.466},
			{153.793, 906.387},
			{158.854, 1045.46},
			{162.331, 1187.25},
			{164.538, 1330.62}
	}; // std devs for hilbert sort

	public SchloeglSystem(double[] c, double[] X0, double tau, double T, double N0) {
		this.c = c;
		this.X0 = X0;
		this.tau = tau;
		this.T = T;
		S = new double[][] { { 1, -1, 1, -1 }, { -1, 1, 0, 0 }, { 0, 0, -1, 1 } };
		init();
		this.N0 = N0;
	}

	public int compareTo(MarkovChainComparable m, int i) {
		if (!(m instanceof SchloeglSystem)) {
			throw new IllegalArgumentException("Can't compare an SchloeglSystem with other types of Markov chains.");
		}
		double mx;

		mx = ((SchloeglSystem) m).X[i];
		return (X[i] > mx ? 1 : (X[i] < mx ? -1 : 0));
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("----------------------------------------------\n");
		sb.append(" SchloeglSystem:\n");
		sb.append("X0 =\t" + "{" + X0[0] + ", " + X0[1] + ", " + (N0 - X0[0] - X0[1]) + "}\n");
		sb.append("c =\t" + "{" + c[0] + ", " + c[1] + ", " + c[2] + "}\n");
		sb.append("T =\t" + T + "\n");
		sb.append("tau =\t" + tau + "\n");
		sb.append("steps =\t" + numSteps + "\n");
		sb.append("----------------------------------------------\n\n");

		return sb.toString();
	}

	@Override
	public double getPerformance() {
		return X[0];
//			return X[1];
	}

	

	@Override
	public double[] getPoint() {
		double[] state01 = new double[N];
		for (int i = 0; i < N; i++)
			state01[i] = getCoordinate(i);
		return state01;
	}

	@Override
	public void computePropensities() {
		double x2 = (N0 - X[0] - X[1]);
		a[0] = 0.5 * c[0] * X[0] * (X[0] - 1.0) * X[1];
		a[1] = c[1] * X[0] * (X[0] - 1.0) * (X[0] - 2.0) / 6.0;
		a[2] = c[2] * x2;
		a[3] = c[3] * X[0];
	}

	public double getCoordinate(int j) {


		return NormalDist.cdf01( (X[j] - means[step][j]) / stdDevs[step][j] );
		}

	


}
