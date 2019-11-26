package umontreal.ssj.markovchainrqmc;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.rng.*;
import umontreal.ssj.util.Chrono;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.hups.*;

/**
 * This class defines a generic Markov chain and provides basic tools to
 * simulate it for a given number of steps or until it stops, and to recover the
 * performance measure. Chains can be cloned, so one can simulate many
 * replicates in parallel.
 * In a concrete subclass, it suffices to implement the three abstract
 * methods  {@link umontreal.ssj.markovchain.MarkovChain.initialState()
 * initialState()},
 * {@link umontreal.ssj.markovchain.MarkovChain.nextStep(RandomStream)
 * nextStep(stream)} and
 * {@link umontreal.ssj.markovchain.MarkovChain.getPerformance()
 * getPerformance()} to get things going.  One would usually 
 * implement subclasses of  @ref MarkovChainComparable or
 * @ref MarkovChainDouble rather than direct subclasses of `MarkovChain`.
 * Some other methods are then needed.
 * It is also *essential* to override the method  `clone`, if the class
 * contains non primitive objects, in order to clone these objects.
 *
 * The methods in this class permit one to simulate the chain over a certain number of steps
 * via Monte Carlo or randomized quasi-Monte Carlo.  Statistics on
 * the performance measure of the chain are computed during these simulations.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class MarkovChain implements Cloneable {

   Chrono timer = Chrono.createForSingleThread();
   public int numSteps;                // Number of steps
   protected boolean stopped = false;

   /**
    * @name Abstract methods
    * @{
    */

   /**
    * Sets the Markov chain to its (deterministic) initial state and
    * initializes the collectors for the performance measure, ready to
    * start a simulation. Must also set variable `stopped` to false.
    */
   public abstract void initialState();

   /**
    * Simulates one more step of the chain, from its current state, using
    * `stream` for the randomness. If the chain stops, this method must
    * set variable `stopped` to true.
    */
   public abstract void nextStep (RandomStream stream);

   /**
    * Returns the performance measure (total or average cost or gain) so
    * far, for the current simulation run.
    */
   public abstract double getPerformance();

   /**
    * @}
    */

   /**
    * @name Other methods
    * @{
    */

   /**
    * Returns a clone of the chain.
    */
   public Object clone() throws CloneNotSupportedException {
        MarkovChain o = null;
        try {
            o = (MarkovChain) super.clone ();
        } catch (CloneNotSupportedException e) {
            System.err.println ("This MarkovChain cannot be cloned");
        }
        return o;
    }

   /**
    * Tells if the chain has stopped.
    */
   public boolean hasStopped() {
        return stopped;
    }

   /**
    * Starts a new simulation run and simulates `numSteps` steps of the
    * Markov chain or until the chain stops, using the given `stream`.
    */
   public void simulSteps (int numSteps, RandomStream stream) {
        initialState ();
        this.numSteps = numSteps;
        int step = 0;
        while (step < numSteps && !hasStopped()){
            nextStep (stream);
            ++step;
        }
    }

   /**
    * Starts a new simulation run and simulates until the stopping time is
    * reached, using the given <tt>stream</tt>. Same as
    * {@link #simulSteps() simulSteps(Integer.MAX_VALUE, stream)}.
    */
   public void simulSteps (RandomStream stream) {
       simulSteps (Integer.MAX_VALUE, stream);
   }

   /**
    * Performs `n` simulation runs of the chain, for `numSteps` steps per
    * run, using the given `stream`. The statistics on the performance for
    * the `n` runs are placed in `statRuns`.
    */
   public void simulRuns (int n, int numSteps, RandomStream stream,
                          Tally statRuns) {
        statRuns.init ();
        for (int i = 0; i < n; i++) {
            simulSteps (numSteps, stream);
            statRuns.add (getPerformance ());
        }
    }

   /**
    * Same as  #simulRuns, except that the stream is first reset to its
    * initial seed and then reset to the first substream at the beginning
    * and to the next substream after each run.
    */
   public void simulRunsWithSubstreams (int n, int numSteps,
                                        RandomStream stream, Tally statRuns) {
        statRuns.init ();
        stream.resetStartStream ();
        for (int i = 0; i < n; i++) {
            simulSteps (numSteps, stream);
            statRuns.add (getPerformance ());
            stream.resetNextSubstream ();
        }
    }

	/**
	 * Perform n simulation runs of the chain, each for numSteps steps, and
	 * returns average.
	 */
	public double simulMC (int n, int numSteps) {
		Tally statRuns = new Tally();
		simulRunsWithSubstreams (n, numSteps, new MRG32k3a(), statRuns);
		return statRuns.average();
	}

	/**
	 * Perform n runs, each one until the chain stops.
	 */
	public double simulMC (int n) {
		  return simulMC (n, Integer.MAX_VALUE);
   }

	/**
	 * Perform n runs and compute the average, reapeat m times and return the
	 * stats in t.
	 */
	public void simulRepMC (int n, int numSteps, int m, Tally t) {
		 for (int rep = 0; rep < m; ++rep) {
 			  t.add (simulMC (n, numSteps));
     } 
	}

	/**
	 * Same as previous one, but run the chains until they stop.
	 */
	public void simulRepMC (int n, int m, Tally t) {
		simulRepMC (n, Integer.MAX_VALUE, m, t);
   }

   /**
    * Performs `m` independent replicates of @f$n@f$ simulation runs of
    * the chain using a RQMC point set, each time storing the average of
    * the performance over the @f$n@f$ chains. @f$n@f$ is the number of
    * points in RQMC point set `p`. Each run goes for `numSteps` steps.
    * For each replicate, the point set `p` is randomized using `rand`, an
    * iterator is created, and each run uses a different substream of this
    * iterator (i.e., a different point). The statistics on the
    * performance for the `m` independent replications are placed in
    * `statReps`.
    */
   public void simulRQMC (PointSet p, int m, int numSteps,
                          PointSetRandomization rand, Tally statReps) {
        statReps.init ();
        Tally statRuns = new Tally ();   // Used within the runs.
        int n = p.getNumPoints();        // Number of points.
        RandomStream stream = p.iterator ();
        for (int rep = 0; rep < m; rep++) {
            p.randomize(rand);
            simulRunsWithSubstreams (n, numSteps, stream, statRuns);
            statReps.add (statRuns.average ());
        }
    }

   /**
    * Same as  #simulRuns but also returns the results as a formatted
    * string.
    */
   public String simulRunsFormat (int n, int numSteps, RandomStream stream,
                                  Tally statRuns) {
        timer.init ();
        simulRuns (n, numSteps, stream, statRuns);
        StringBuffer sb = new StringBuffer
           ("----------------------------------------------" +
                PrintfFormat.NEWLINE);
        sb.append ("MC simulations:" + PrintfFormat.NEWLINE);
        sb.append (" Number of runs n          = " + n  +
                PrintfFormat.NEWLINE);
        sb.append (formatResults (statRuns));
        sb.append (" CPU Time = " + timer.format () + PrintfFormat.NEWLINE);
        return sb.toString ();
    }

   /**
    * Same as  #simulRunsWithSubstreams but also returns the results as a
    * formatted string.
    */
   public String simulRunsWithSubstreamsFormat (int n, int numSteps,
                                                RandomStream stream,
                                                Tally statRuns) {
        timer.init ();
        simulRunsWithSubstreams (n, numSteps, stream, statRuns);
        StringBuffer sb = new StringBuffer
           ("----------------------------------------------" +
             PrintfFormat.NEWLINE);
        sb.append ("MC simulations with substreams:" + PrintfFormat.NEWLINE);
        sb.append (" Number of runs n          = " + n  +
           PrintfFormat.NEWLINE);
        sb.append (formatResults (statRuns));
        sb.append (" CPU Time = " + timer.format () + PrintfFormat.NEWLINE);
        return sb.toString ();
    }

   /**
    * Same as  #simulRQMC but also returns the results as a formatted
    * string.
    */
   public String simulRQMCFormat (PointSet p, int m, int numSteps,
                                  PointSetRandomization rand, Tally statReps) {
        timer.init();
        simulRQMC (p, m, numSteps, rand, statReps);
        int n = p.getNumPoints();
        StringBuffer sb = new StringBuffer
            ("----------------------------------------------" +
                PrintfFormat.NEWLINE);
        sb.append ("RQMC simulations:" + PrintfFormat.NEWLINE +
                PrintfFormat.NEWLINE);
        sb.append (p.toString ());
        sb.append (PrintfFormat.NEWLINE + " Number of indep. randomization, m = "
             + m  + PrintfFormat.NEWLINE);
        sb.append (" Number of points n        = "+ n  +
                    PrintfFormat.NEWLINE);
        sb.append (formatResultsRQMC (statReps, n));
        sb.append (" CPU Time = " + timer.format () + PrintfFormat.NEWLINE);
        return sb.toString ();
    }

   /**
    * Similar to  #simulRQMCFormat, but also gives the variance
    * improvement factor with respect to MC. Assuming that `varMC` gives
    * the variance per run for MC.
    */
   public String testImprovementRQMCFormat (PointSet p, int m, int numSteps,
                                      PointSetRandomization rand, double varMC,
                                      Tally statReps) {
      // Removed next line because numSteps may be infinite!
      // p.randomize (0, numSteps * dimPerStep, noise);
      StringBuffer sb = new StringBuffer (simulRQMCFormat
              (p, m, numSteps, rand, statReps));
      double var = p.getNumPoints() * statReps.variance();
      sb.append (" Variance ratio: " +
           PrintfFormat.format (15, 10, 4, varMC/var) +
                PrintfFormat.NEWLINE);
      return sb.toString ();
}

   /**
    * Returns a string containing the mean, the variance, and a 90%
    * confidence interval for `stat`.
    */
   public String formatResults (Tally stat) {
        StringBuffer sb = new StringBuffer (" Average value             = ");
        sb.append (PrintfFormat.format (12, 9, 5, stat.average ()) +
                   PrintfFormat.NEWLINE);
        sb.append (" Variance                  = ");
        sb.append (PrintfFormat.format (12, 9, 5, stat.variance ()) +
                   PrintfFormat.NEWLINE);
        sb.append (stat.formatCIStudent (0.9, 7));
        return sb.toString ();
    }

   /**
    * Returns a string containing the mean, the variance multiplied by
    * `numPoints`, and a 90% confidence interval for `stat`.
    */
   public String formatResultsRQMC (Tally stat, int numPoints) {
        StringBuffer sb = new StringBuffer (" Average value             = ");
        sb.append (PrintfFormat.format (12, 9, 5, stat.average ()) +
                   PrintfFormat.NEWLINE);
        sb.append (" Variance * numPoints      = ");
        sb.append (PrintfFormat.format (12, 9, 5, numPoints * stat.variance ()) +
                   PrintfFormat.NEWLINE);
        sb.append (stat.formatCIStudent (0.9, 7));
        return sb.toString ();
    }

}

/**
 * @}
 */