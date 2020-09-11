package umontreal.ssj.markovchainrqmc.florian.chemicalReactionNetworks;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import umontreal.ssj.markovchainrqmc.MarkovChainComparable;
import umontreal.ssj.probdist.PoissonDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.sort.florian.MultiDim;

/**
 * Implements an abstract class to simulate chemical reaction networks with
 * Gillespie's @f$\tau@f$ algorithm. These aim to observe the number of
 * molecules of @f$\ell@f$ chemical species @f$S_1,S_2,\dots,S_{\ell}@f$ that
 * can react via @f$K@f$ predefined chemical reactions @f$R_1,R_2,\dots,R_K@f$
 * over a time interval @f$[0,T]@f$. This is basically done by simulating a
 * discrete time Markov chain. More precisely, one subdivides the time interval
 * into @f$d@f$ smaller intervals of length @f$au>0@f$. Let us denote the copy
 * numbers at time @fj\tau@f$, i.e., after the @f$j@f$th step by $\mathbf X_j$.
 * Then, the system is updated via
 * 
 * @f[ \mathbf{X}_{j+1} = \mathbf{X}_j + \cdot \mathbf{p}, @f]
 * 
 * where @f$\mathbf{S}@f$ denotes the @f$\ell\times K@f$ stoichiometric matrix,
 * whose entry at position @f$(i,k)@f$ gives how many molecules of
 * type @f$S_i@f$ are created/lost by one reaction of type @f$R_k@f$; and where
 * the vector @f$\mathbf{p}=(p_1,p_2,\dots,p_K)@f$ contains @f$K@f$ poisson
 * variates simulating how often each reaction fires within the time
 * interval @f$[j\tau,(j+1)\tau)@f$.
 * 
 * The mean parameters of the poisson variate @f$p_k@f$ are given by
 * $a_k(\mathbf{X}_j)\tau$, where @f$a_k@f$ gives the propensity function of
 * reaction @f$R_k@f$. Since the actual choice of propensity functions involve a
 * modeling choice and is different for each example, they are implemented as
 * abstract functions.
 * 
 * @author florian
 *
 */
public abstract class ChemicalReactionNetwork extends MarkovChainComparable implements MultiDim {

	/**
	 * Current step.
	 */
	int step;
	/**
	 * Initial data
	 */
	double[] X0;
	/**
	 * Propensity functions
	 */
	double[] a;
	/**
	 * Stoichiometric matrix.
	 */
	double[][] S;
	/**
	 * Reaction rates.
	 */
	double[] c;

	/**
	 * Number of reactions in the system.
	 */
	public int K;

	/**
	 * Time step.
	 */
	double tau;
	/**
	 * Final time.
	 */
	double T;
	/**
	 * Copy numbers, i.e., the states of the chain.
	 */
	double[] X;
	/**
	 * Number of molecular species.
	 */
	int N;

	/**
	 * Initializes/resets several parameters to their initial states.
	 */
	public void init() {
		K = c.length;
		N = X0.length;
		X = new double[N];
		a = new double[K];
		numSteps = (int) Math.ceil(T / tau);
		stateDim = N;
		X = new double[N];
		for (int j = 0; j < N; ++j)
			X[j] = X0[j];
//		initialState();
	}

	/**
	 * Setter for X0.
	 * 
	 * @param X0
	 */
	public void setInitialState(double[] X0) {
		this.X0 = X0;
	}

	public int getK() {
		return K;
	}

	// Initial value of X.
	public void initialState() {
		step = 0;
		for (int j = 0; j < N; ++j)
			X[j] = X0[j];
	}

	/**
	 * Propensity functions.
	 */
	abstract public void computePropensities();

	/**
	 * Simulates the next step \a reps times and takes the the average over these
	 * repetitions as the next state.
	 * 
	 * @remark **Florian:** If we want to use this with RQMC points, we need to find
	 *         a way to reset the stream to the previous point.
	 * 
	 * @param stream RandomStream used.
	 * @param state the current state.  
	 * @param reps the number of repetitions.
	 */
	public void nextStep(RandomStream stream, double[] state, int reps) {
		// write state before step
		double[] tmp = getPoint();
		for (int d = 0; d < tmp.length; d++)
			state[d] = tmp[d];

		double[] futureState = new double[state.length];
		nextStep(stream);
		for (int j = 0; j < state.length; j++)
			futureState[j] = (getPoint())[j];

		for (int r = 1; r < reps; r++) {
			// TODO: for sobol, etc. need to reset stream to prev. point.

			// reset X to old state
			for (int j = 0; j < state.length; j++)
				X[j] = state[j];

			// carry out next step
			step--;
			nextStep(stream);

			// update the future state avg.
			for (int j = 0; j < state.length; j++) {
				futureState[j] = ((double) r) * futureState[j] + (getPoint())[j];
				futureState[j] /= (double) (r + 1.0);
//				futureState[j] += (getState())[j];
			}

		}
		// after all the repetitions, set the state to the avg
		for (int j = 0; j < state.length; j++)
			X[j] = futureState[j] ;// / (double) reps;
	}

	/**
	 * Simulate \a numSteps steps of the chain with using the #nextStep method that
	 * repeats each step \a reps times and takes the average of these repetitions as
	 * the future state.
	 * 
	 * @param numSteps number of steps to be simulated
	 * @param stream   the random stream used.
	 * @param states   the state at each step
	 * @param reps     the number of repetitions
	 */
	public void simulSteps(int numSteps, RandomStream stream, double[][] states, int reps) {
		initialState();
		this.numSteps = numSteps;
		int step = 0;
		while (step < numSteps && !hasStopped()) {
			states[step] = new double[getPoint().length];
			nextStep(stream, states[step], reps);
			++step;
		}
	}

	/**
	 * Simulates \a n chains over \a numSteps steps using the #nextStep method that
	 * repeats each step \a reps times and takes the average of these repetitions as
	 * the future state.
	 * @param n number of chains.
	 * @param numSteps number of steps
	 * @param stream the random stream used.
	 * @param states the states of each chain at each step.
	 * @param performance array to which the performance of each chain is written.
	 * @param reps the number of replications.
	 */
	public void simulRuns(int n, int numSteps, RandomStream stream, double[][][] states, double[] performance,
			int reps) {
		for (int i = 0; i < n; i++) {
			states[i] = new double[numSteps][];
			simulSteps(numSteps, stream, states[i], reps);
			performance[i] = getPerformance();
		}
	}



	@Override
	public void nextStep(RandomStream stream) {
		step++;
		double[] p = new double[K];
		computePropensities();


		for (int k = 0; k < K; k++) {
			p[k] = PoissonDist.inverseF(a[k] * tau, stream.nextDouble());
		}

		double[] temp = new double[N];
		for (int j = 0; j < N; j++)
			temp[j] = X[j];
		for (int n = 0; n < N; n++) {

			for (int k = 0; k < K; k++) {
				temp[n] += S[n][k] * p[k];
			}


		}
		X = temp;


	}

	@Override
	public int dimension() {
		return N;
	}

	
	public double[] getPoint() {
		return X;
	}

	/**
	 * Simulates \a n chains over \a numSteps steps using the #nextStep method that
	 * repeats each step \a reps times and takes the average of these repetitions as
	 * the future state. The data is displayed and written to a csv-file at the location \a dataPath. 
	 * @param dataPath location of the file.
	 * @param dataLabel short description as part of the file name.
	 * @param n number of chains.
	 * @param numSteps number of steps.
	 * @param stream the random stream used.
	 * @param reps the number of repetitions.
	 * @throws IOException
	 */
	public void genData(String dataPath, String dataLabel, int n, int numSteps, RandomStream stream, int reps)
			throws IOException {
		double[][][] states = new double[n][][];
		double[] performance = new double[n];
		simulRuns(n, numSteps, stream, states, performance, reps);
		StringBuffer sb;
		FileWriter fw;
		File file;
		for (int step = 0; step < numSteps; step++) {
			sb = new StringBuffer("");
			file = new File(dataPath + dataLabel + "_Step_" + step + ".csv");
//			file.getParentFile().mkdirs();
			fw = new FileWriter(file);

			for (int i = 0; i < n; i++) {
				for (int j = 0; j < getStateDimension(); j++)
					sb.append(states[i][step][j] + ",");
				sb.append(performance[i] + "\n");
			}
			fw.write(sb.toString());
			fw.close();
			System.out.println("*******************************************");
			System.out.println(" STEP " + step);
			System.out.println("*******************************************");
			System.out.println(sb.toString());
		}

	}

	
}
