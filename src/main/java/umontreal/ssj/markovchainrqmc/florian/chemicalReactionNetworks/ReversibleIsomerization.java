package umontreal.ssj.markovchainrqmc.florian.chemicalReactionNetworks;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.util.Chrono;

import umontreal.ssj.util.sort.MultiDim01;
import umontreal.ssj.util.sort.florian.MultiDim;

/**
 * Implements the Reversible Isomerization example from, e.g., Beentjes and
 * Baker '18. Since the total number of molecules does not change through any
 * reaction, the state space is, in fact, one dimensional.
 * 
 * @author florian
 *
 */
public class ReversibleIsomerization extends ChemicalReactionNetwork implements MultiDim {

	/**
	 * Total number of molecules in the system
	 */
	double N0; 

	/**
	 * Constructor
	 * @param c reaction rates.
	 * @param X0 initial state.
	 * @param tau step size.
	 * @param T final time.
	 * @param N0 Total number of molecules in the system.
	 */
	public ReversibleIsomerization(double[] c, double[] X0, double tau, double T, double N0) {
		this.c = c;
		this.X0 = X0;
		this.tau = tau;
		this.T = T;
		S = new double[][] { { -1, 1 }, { 1, -1 } };
		init();
		this.N0 = N0;
	}

	public int compareTo(MarkovChainComparable m, int i) {
		if (!(m instanceof ReversibleIsomerization)) {
			throw new IllegalArgumentException("Can't compare an ReversibleIso with other types of Markov chains.");
		}
		double mx;

		mx = ((ReversibleIsomerization) m).X[i];
		return (X[i] > mx ? 1 : (X[i] < mx ? -1 : 0));
	}

	public String toString() {
		StringBuffer sb = new StringBuffer("----------------------------------------------\n");
		sb.append(" ReversibleIsomerisation:\n");
		sb.append("X0 =\t" + "{" + X0[0] + ", " + (N0 - X0[0]) + "}\n");
		sb.append("c =\t" + "{" + c[0] + ", " + c[1] + "}\n");
		sb.append("T =\t" + T + "\n");
		sb.append("tau =\t" + tau + "\n");
		sb.append("steps =\t" + numSteps + "\n");
		sb.append("----------------------------------------------\n\n");

		return sb.toString();
	}

	@Override
	public void computePropensities() {
		double x1 = (N0 - X[0]);
		a[0] = c[0] * X[0];
		a[1] = c[1] * x1;
	}

	@Override
	public double getPerformance() {
		return X[0];
	}

	@Override
	public double[] getPoint() {
		double[] state01 = new double[N];
		for (int i = 0; i < N; i++)
			state01[i] = getCoordinate(i);
		return state01;
	}

	public double getCoordinate(int j) {
		return X[0];
	}

	


}
