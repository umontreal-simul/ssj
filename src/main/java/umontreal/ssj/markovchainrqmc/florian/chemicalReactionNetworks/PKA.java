package umontreal.ssj.markovchainrqmc.florian.chemicalReactionNetworks;

import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.util.sort.MultiDim01;

/**
 * Implements the #ChemicalReactionNetwork for the cyclic AMP activation of PKA,
 * see Koh and Blackwell '12, Strehl and Ilie '15.
 * 
 * @author florian
 *
 */
public class PKA extends ChemicalReactionNetwork implements MultiDim01 {

	/**
	 * Constructor
	 * @param c reaction rates
	 * @param X0 initial states
	 * @param tau step length
	 * @param T final time
	 */
	public PKA(double[] c, double[] X0, double tau, double T) {
		this.c = c;
		this.X0 = X0;
		this.tau = tau;
		this.T = T;
		S = new double[][] { { -1, 1, 0, 0, 0, 0 }, { -2, 2, -2, 2, 0, 0 }, { 1, -1, -1, 1, 0, 0 },
				{ 0, 0, 1, -1, -1, 1 }, { 0, 0, 0, 0, 1, -1 }, { 0, 0, 0, 0, 2, -2 } };
		init();
	}

	/**
	 * Means of the states over 15 steps for the example in our paper
	 */
	double[][] means = { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
			{ 27782.2, 22132.8, 6087.02, 1334.56, 1096.23, 1092.47 },
			{ 25809.7, 17041., 7486.09, 1911.67, 1092.53, 1085.06 },
			{ 24723.3, 14032.1, 8154.43, 2333.36, 1088.89, 1077.78 },
			{ 24017.7, 12003.3, 8551.22, 2645.74, 1085.31, 1070.62 },
			{ 23516.2, 10526.3, 8815.82, 2886.21, 1081.79, 1063.58 },
			{ 23138.5, 9395.29, 9005.64, 3077.53, 1078.32, 1056.65 },
			{ 22842.5, 8497.52, 9148.81, 3233.8, 1074.92, 1049.84 },
			{ 22603.4, 7765.33, 9260.83, 3364.19, 1071.56, 1043.12 },
			{ 22405.9, 7155.42, 9350.99, 3474.89, 1068.26, 1036.51 },
			{ 22239.6, 6638.78, 9425.19, 3570.2, 1065.01, 1030.02 },
			{ 22097.5, 6194.86, 9487.37, 3653.29, 1061.81, 1023.62 },
			{ 21974.6, 5809.01, 9540.26, 3726.46, 1058.66, 1017.32 },
			{ 21867.1, 5470.24, 9585.82, 3791.48, 1055.55, 1011.11 },
			{ 21772.3, 5170.25, 9625.49, 3849.69, 1052.5, 1005. } };

	/**
	 * Standard deviations of the states over 15 steps
	 */
	double[][] stdDevs = { { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, { 72.2169, 147.578, 73.8083, 15.3215, 1.94182, 3.88363 },
			{ 60.9309, 130.231, 68.0601, 27.0314, 2.71369, 5.42738 },
			{ 54.8053, 114.932, 67.8547, 30.9913, 3.27786, 6.55572 },
			{ 51.2789, 103.649, 68.6454, 32.9158, 3.7385, 7.47699 },
			{ 49.003, 95.0904, 69.5041, 34.0796, 4.13071, 8.26142 },
			{ 47.4479, 88.453, 70.2515, 34.8499, 4.467, 8.934 },
			{ 46.3324, 83.1412, 70.8924, 35.4005, 4.76759, 9.53519 },
			{ 45.4829, 78.7204, 71.4452, 35.8304, 5.03575, 10.0715 },
			{ 44.8055, 74.932, 71.8831, 36.1387, 5.27811, 10.5562 },
			{ 44.2434, 71.69, 72.2576, 36.4183, 5.49952, 10.999 },
			{ 43.7931, 68.8705, 72.5846, 36.6378, 5.70143, 11.4029 },
			{ 43.4126, 66.3247, 72.8685, 36.8135, 5.8875, 11.775 },
			{ 43.0824, 64.06, 73.1321, 36.9902, 6.06264, 12.1253 },
			{ 42.7897, 62.0298, 73.3465, 37.1401, 6.22393, 12.4479 } };

	
	public String toString() {
		StringBuffer sb = new StringBuffer("----------------------------------------------\n");
		sb.append(" cAMP activation of PKA:\n");
		sb.append("Number of reactions K = " + K + "\n");
		sb.append("Number of species N = " + N + "\n");
		sb.append("X0 =\t" + "{" + X0[0]);
		for (int i = 1; i < X0.length; i++)
			sb.append(", " + X0[i]);
		sb.append("}\n");

		sb.append("c =\t" + "{" + c[0]);
		for (int i = 1; i < c.length; i++)
			sb.append(", " + c[i]);
		sb.append("}\n");
		sb.append("T =\t" + T + "\n");
		sb.append("tau =\t" + tau + "\n");
		sb.append("steps =\t" + numSteps + "\n");
		sb.append("----------------------------------------------\n\n");

		return sb.toString();
	}

	@Override
	public int compareTo(MarkovChainComparable m, int i) {
		if (!(m instanceof PKA)) {
			throw new IllegalArgumentException("Can't compare a PKA with other types of Markov chains.");
		}
		double mx;

		mx = ((PKA) m).X[i];
		return (X[i] > mx ? 1 : (X[i] < mx ? -1 : 0));
	}

	@Override
	public double[] getPoint() {
		double[] state01 = new double[N];
		for (int i = 0; i < N; i++)
			state01[i] = getCoordinate(i);
		return state01;
	}

	@Override
	public double getCoordinate(int j) {
		double zvalue;


		return NormalDist.cdf01((X[j] - means[step][j]) / stdDevs[step][j]);
	}

	@Override
	public void computePropensities() {
		a[0] = 0.5 * c[0] * X[0] * X[1] * (X[1] - 1.0);
		a[1] = c[1] * X[2];
		a[2] = 0.5 * c[2] * X[2] * X[1] * (X[1] - 1.0);
		a[3] = c[3] * X[3];
		a[4] = c[4] * X[3];
		a[5] = 0.5 * c[5] * X[4] * X[5] * (X[5] - 1.0);
	}

	@Override
	public double getPerformance() {
		return X[0]; // PKA
//		return X[1]; //cAMP
//		return X[2]; //PKA-cAMP2
//		return X[3];
//		return X[4]; //PKAr
//		return X[5]; //PKAc
	}

}
