package umontreal.ssj.markovchainrqmc;
 import umontreal.ssj.stat.Tally;
 import umontreal.ssj.rng.RandomStream;
 import umontreal.ssj.util.*;
 import umontreal.ssj.util.sort.*;
 import umontreal.ssj.hups.*;
 import umontreal.ssj.charts.*;
 import umontreal.ssj.functionfit.LeastSquares;
 import java.util.*;
 import java.io.IOException;
 import java.io.Writer;
 import java.io.FileWriter;
 import java.lang.reflect.Array;

/**
 * This class provides tools to simulate an array of
 * @ref MarkovChainComparable objects with the array-RQMC method of
 * @cite vLEC08a, @cite vLEC09d&thinsp;. It offers tools to construct the
 * array of chains, simulate with array-RQMC for one step or multiple steps,
 * perform experiments, and report some results.
 *
 * At each step, the @f$n@f$ states (which correspond to the @f$n@f$ copies
 * of the chain) are sorted in @f$\ell@f$ dimensions (usually the dimension
 * of the state). If @f$\ell=1@f$, this is an ordinary sort. If @f$\ell>
 * 1@f$, they are sorted using a  @ref umontreal.ssj.util.MultiDimSort, which
 * in general can sort  @ref umontreal.ssj.util.MultiDimComparable objects.
 * As a special case, they can also be sorted using a
 * @ref umontreal.ssj.util.MultiDimSort01 type of sort. Some of these sorts
 * effectively map the chain states to the one-dimensional interval
 * @f$(0,1)@f$. For example, the sorts based on a Hilbert curve, such as
 * @ref umontreal.ssj.util.HilbertCurveSort, do that. In that case, we put
 * @f$\ell’=1@f$, otherwise @f$\ell’=\ell@f$.
 *
 * To move the chains ahead, we use an @f$(\ell’+d)@f$-dimensional RQMC
 * @ref umontreal.ssj.hups.PointSet. The first @f$\ell’@f$ coordinates are
 * used to sort the points to map them to the (sorted) chains. In some cases,
 * when @f$\ell’=1@f$, the first coordinate of the points is not stored
 * explicitly; the points are already sorted by a first coordinate which is
 * only implicit (it serves to enumerate the points). This is what happens
 * for a Sobol’ sequence, or a  @ref umontreal.ssj.hups.KorobovLattice for
 * which the first coordinate is dropped, for example. Only @f$d@f$
 * coordinates are produced with those point sets and have to be randomized.
 * In the other cases, the points must be sorted by their first @f$\ell’@f$
 * coordinates.
 *
 * At each step the  @ref umontreal.ssj.hups.PointSet is randomized by using
 * a  @ref umontreal.ssj.hups.PointSetRandomization. There are types of point
 * sets for which only the last @f$d@f$ coordinates have to be randomized
 * (e.g., digital nets and lattice rules) and others for which all @f$\ell’
 * + d@f$ coordinates must be randomized (e.g., a stratified sample).
 *
 * When applying the array-RQMC method, the number of point coordinates used
 * explicitly to sort the points and match the chains (which is either
 * @f$\ell’@f$ or 0) is passed in a variable named `sortCoord`, and the
 * number of these coordinates that must be randomized at each step is passed
 * in a variable named `sortCoordRand`. In the case where @f$\ell’=1@f$ and
 * the first coordinate is only implicit, both should be 0. In the case where
 * @f$\ell’ > 1@f$ but only the last @f$d@f$ coordinates have to be
 * randomized at each step, one would put `sortCoord` to @f$\ell’@f$ and
 * `sortCoordRand` to 0.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class ArrayOfComparableChains <T extends MarkovChainComparable> {
   protected T baseChain;  // The base chain.
   protected int n;               // Current number of chains.
   // protected int stateDim;        // Dimension of the chain state.
   protected T[] chains;  // Array of n comparable chains.
   protected double[] performances;  // Performances for the n chains.
   protected PointSetRandomization randomization;
   protected MultiDimSort<T> savedSort;
	protected int sortCoordPts = 0;   // Point coordinates used to sort points.

   /**
    * Creates an array of the comparable chain `baseChain`. The method
    * #makeCopies(int) must be called afterward to make the actual copies
    * of the chain.
    */
   public ArrayOfComparableChains (T baseChain) {
      this.baseChain = baseChain;
      // stateDim = baseChain.stateDim;
   }

   /**
    * Creates an array of the comparable chain `baseChain`. The method
    * #makeCopies(int) must be called to make the copies. `rand` will be
    * used to randomize the point sets in the simulations. `sort` will be
    * used to sort the chains.
    */
   public ArrayOfComparableChains (T baseChain,
                                   PointSetRandomization rand,
                                   MultiDimSort sort) {
      this.baseChain = baseChain;
      // stateDim = baseChain.stateDim;
      randomization = rand;
      savedSort = sort;
   }

   /**
    * Creates <tt>n</tt> copies (clones) of the chain <tt>baseChain</tt>
    * and puts them in an array, ready for the array RQMC simulation.
    */
   public void makeCopies (int n) {

		@SuppressWarnings("unchecked")
		final T[] c = (T[]) Array.newInstance(baseChain.getClass(), n);
		chains = c; // Array of Markov chains.

      this.n = n;
      performances = new double[n];  // Array to store the performances.
      for (int i = 0; i < n; i++) {
         try {
				chains[i] = (T) baseChain.clone ();
         }
			catch (CloneNotSupportedException e) {
            System.err.println ("ArrayOfComparableChains:");
            e.printStackTrace ();
         }
      }
   }

   /**
    * Initializes the <tt>n</tt> copies (clones) of the chain
    * <tt>baseChain</tt> to their initial state by calling
    * {@link umontreal.ssj.markovchain.MarkovChain.initialState()
    * initialState()} on each chain.
    */
   public void initialStates() {
      int i = 0;
      for (T mc : chains){
         mc.initialState();
         performances[i] = mc.getPerformance();  // Needed?   Why do this? 
         ++i;
       }
   }

   /**
    * Returns the number `n` of chains.
    */
   public int getN() {
       return n;
   }

   /**
    * Returns the underlying array of `n`  @ref MarkovChainComparable.
    */
   public T[] getChains() {
       return chains;
   }

   /**
    * Sets the internal  @ref umontreal.ssj.hups.PointSetRandomization to
    * `rand`.
    */
   public void setRandomization (PointSetRandomization rand) {
       randomization = rand;
   }

   /**
    * Returns the internal  @ref umontreal.ssj.hups.PointSetRandomization.
    */
   public PointSetRandomization getRandomization() {
       return randomization;
   }

   /**
    * Sets the internal  @ref umontreal.ssj.util.MultiDimSort to `sort`.
    */
   public void setSort (MultiDimSort sort) {
       savedSort = sort;
   }

   /**
    * Returns the saved  @ref umontreal.ssj.util.MultiDimSort.
    */
   public MultiDimSort getSort () {
       return savedSort;
   }

   /**
    * Randomized the point set `p` and Simulates the @f$n@f$ copies of the
    * chain, one step for each copy, using
    * @ref umontreal.ssj.hups.PointSet `p`, where @f$n@f$ is the current
    * number of copies (clones) of the chain and is *assumed* to equal the
    * number of points in `p`. The points are randomized before the
    * simulation using the stored
    * @ref umontreal.ssj.hups.PointSetRandomization. If `sortCoordPts`
    * @f$>0@f$, the points are also sorted explicitly based on their first
    * `sortCoordPts` coordinates, at each step, after they are randomized.
    * In that case, `p` must implement the
    * @ref umontreal.ssj.util.MultiDim01 interface. The dimension of `p`
    * must be at least as large as `sortCoordPts` @f$+@f$ the number of
    * uniforms required to simulate one step of the chain. Returns the
    * number of chains that have not stopped yet.
    */
   public int simulOneStepArrayRQMC (PointSet p, PointSetRandomization rand, 
	       MultiDimSort sort, int sortCoordPts) {
      int nStopped = 0;
      p.randomize(rand);                        // Randomize point set.
         if (sortCoordPts > 0) { 
             if (!(p instanceof CachedPointSet))
                throw new IllegalArgumentException("p is not a CachedPointSet.");
				 if (sortCoordPts > 1)
                ((CachedPointSet) p).sort(sort);   // Sort points using first sortCoordPts coordinates. 
             else
                ((CachedPointSet) p).sortByCoordinate (0);  // Sort by first coordinate.
         }
      PointSetIterator stream = p.iterator ();
      stream.resetCurPointIndex ();             // Go to first point.
      int i = 0;
      for (T mc : chains) { // Assume the chains are sorted
         if (mc.hasStopped()) {
            ++nStopped;
         } else {
            stream.setCurCoordIndex (sortCoordPts); // Skip first sortCoordPts coord.
            mc.nextStep (stream);               // simulate next step of the chain.
            stream.resetNextSubstream ();       // Go to next point.
            if (mc.hasStopped())
               ++nStopped;
         }
         performances[i] = mc.getPerformance();
         ++i;
      }
      return n - nStopped;
    }

   /**
    * This version uses the preselected randomization and sort, with
    * `sortCoordPts = 0`.
    */
   public int simulOneStepArrayRQMC (PointSet p) {
      return simulOneStepArrayRQMC (p, randomization, savedSort, 0);
   }

   /**
    * Simulates the @f$n@f$ copies of the chain, `numSteps` steps for each
    * copy, using  @ref umontreal.ssj.hups.PointSet `p`, where @f$n@f$ is
    * the current number of copies (clones) of the chain and is *assumed*
    * to equal the number of points in `p`. At each step, the points are
    * randomized using `rand`. All coordinates are randomized. If
    * `sortCoordPts` @f$>0@f$, the points are also sorted explicitly based
    * on their first `sortCoordPts` coordinates, at each step, after they
    * are randomized. In that case, `p` must implement the
    * @ref umontreal.ssj.util.MultiDim01 interface. If the coordinates
    * used for the sort do not have to be randomized at each step and the
    * points do not have to be sorted again, one should remove these
    * coordinates before invoking this method and use `sortCoordPts=0`. In
    * this case, the points must be sorted before invoking this method.
    * The class  @ref umontreal.ssj.hups.SortedAndCutPointSet can be
    * useful for this. The dimension of `p` must be at least as large as
    * `sortCoordPts` @f$+@f$ the number of uniforms required to simulate
    * one step of the chain. The method returns the average performance
    * per run. An array that contains the performance for each run can
    * also be obtained via  #getPerformances()(.)
    */
   public double simulArrayRQMC (PointSet p, PointSetRandomization rand, 
	       MultiDimSort sort, int sortCoordPts, int numSteps) {
      int numNotStopped = n;
      initialStates();
      int step = 0;
      while (step < numSteps && numNotStopped > 0) {
         if (numNotStopped == n) sort.sort(chains, 0, n); 
				 else sortNotStoppedChains (sort); // Sort the numNotStopped first chains.          
         p.randomize(rand);           // Randomize the point set.
         if (sortCoordPts > 0) { 
             if (!(p instanceof CachedPointSet))
                throw new IllegalArgumentException("p is not a CachedPointSet.");
				 if (sortCoordPts > 1)
                ((CachedPointSet) p).sort(sort);   // Sort points using first sortCoordPts coordinates. 
             else
                ((CachedPointSet) p).sortByCoordinate (0);  // Sort by first coordinate.
         }
         PointSetIterator stream = p.iterator ();
         stream.resetCurPointIndex ();             // Go to first point.
         int i = 0;
         for (T mc : chains) { // Assume the chains are sorted
            if (mc.hasStopped()) {
               numNotStopped--;
            } else {
	       stream.setCurCoordIndex (sortCoordPts); // Skip first sortCoordPts coord.
               mc.nextStep (stream);                // simulate next step of the chain.
               stream.resetNextSubstream ();        // Go to next point.
               if (mc.hasStopped())
                  numNotStopped--;
            }
            performances[i] = mc.getPerformance();
            ++i;
         }
         ++step;
      }
      return calcMeanPerf();
   }

   /**
    * This version assumes that `sortCoordPts = 0`, so that there is no
    * need to sort the points at each step.
    */
   public double simulArrayRQMC (PointSet p, PointSetRandomization rand, 
         MultiDimSort sort, int numSteps) {
      return simulArrayRQMC (p, rand, sort, 0, numSteps);
   }

   /**
    * This version assumes that `sortCoordPts = 0` and uses the preset
    * randomization and sort for the chains.
    */
   public double simulArrayRQMC (PointSet p, int numSteps) {
      return simulArrayRQMC (p, randomization, savedSort, 0, numSteps);
   }

   /**
    * Returns the vector for performances for the @f$n@f$ chains.
    */
   public double[] getPerformances() {
        return performances;
    }

   /**
    * Computes and returns the mean performance of the @f$n@f$ chains.
    */
   public double calcMeanPerf() {
        double sumPerf = 0.0;                     // Sum of performances.
        for (int i=0; i<n; ++i) {
            sumPerf += performances[i];
        }
        return sumPerf/n;
    }

   /**
    * Performs <tt>m</tt> independent replications of an array-RQMC
    * simulation as in <tt>simulArrayRQMC</tt>. The statistics on the
    * <tt>m</tt> corresponding averages are collected in
    * <tt>statReps</tt>.
    */
   public void simulReplicatesArrayRQMC (PointSet p, PointSetRandomization rand, 
	        MultiDimSort sort, int sortCoordPts, 
          int numSteps, int m, Tally statReps) {
      makeCopies (p.getNumPoints());
      statReps.init ();
      for (int rep = 0; rep < m; rep++) {
         statReps.add (simulArrayRQMC (p, rand, sort, sortCoordPts, numSteps));
      }
   }

   /**
    * Performs <tt>m</tt> independent replications of an array-RQMC
    * simulation as in <tt>simulFormatArrayRQMC</tt>. The statistics on
    * the <tt>m</tt> corresponding averages are collected in
    * <tt>statReps</tt> and the results are also returned in a string.
    */
   public String simulReplicatesArrayRQMCFormat (PointSet p, PointSetRandomization rand, 
	        MultiDimSort sort, int sortCoordPts, 
          int numSteps, int m, Tally statReps) {
      Chrono timer = Chrono.createForSingleThread();
      makeCopies (p.getNumPoints());
      timer.init ();
      statReps.init ();
      for (int rep = 0; rep < m; rep++) {
         statReps.add (simulArrayRQMC (p, rand, sort, sortCoordPts, numSteps));
      }
      StringBuffer sb = new StringBuffer
          ("----------------------------------------------" +
            PrintfFormat.NEWLINE);
      sb.append ("Array-RQMC simulations:" + PrintfFormat.NEWLINE);
      sb.append (PrintfFormat.NEWLINE + p.toString() + ":" +
                 PrintfFormat.NEWLINE);
      sb.append (" Number of indep copies m  = " + m);
      sb.append (PrintfFormat.NEWLINE + " Number of points n        = "
          + n  + PrintfFormat.NEWLINE);
      sb.append (baseChain.formatResultsRQMC (statReps, n));
      sb.append (" CPU Time = " + timer.format () + PrintfFormat.NEWLINE);
      return sb.toString ();
    }

   /**
    * Returns a string that reports the the ratio of MC variance per run
    * `varMC` over the RQMC variance per run `varRQMC`.
    */
   public String varianceImprovementFormat (double varRQMC, double varMC) {
      // double varRQMC = p.getNumPoints() * statReps.variance();
      StringBuffer sb = new StringBuffer (" Variance ratio MC / RQMC: " +
           PrintfFormat.format (15, 10, 4, varMC/varRQMC)   + PrintfFormat.NEWLINE);
      return sb.toString ();
   }

	/**
	 * Performs an experiment to estimate the convergence rate of the RQMC
	 * variance as a function of @f$n@f$, by invoking
	 * #simulReplicatesArrayRQMC(r) epeatedly with a given array of point sets
	 * of different sizes @f$n@f$. Returns a string that reports the mean,
	 * variance, and variance reduction factor (VRF) with respect to MC for
	 * each point set, and the estimated convergence rate of the RQMC variance
	 * as a function of @f$n@f$. Assumes that <tt>varMC</tt> is the variance
	 * per run for MC. The string `methodLabel` should be a brief descriptor of
	 * the method (e.g., the type of point set, type of randomization, and type
	 * of sort). If the string `filenamePlot != null`, then the method also
	 * creates a `.tex` file with that name that contains a plot in log scale
	 * of the variance vs @f$n@f$.
	 */
	public String testVarianceRateFormat (PointSet[] pointSets, PointSetRandomization rand, 
	        MultiDimSort sort, int sortCoordPts, 
          int numSteps, int m, double varMC, String filenamePlot, String methodLabel) {

	  int numSets = pointSets.length; // Number of point sets.
		Tally statPerf = new Tally ("Performance");
		double[] logn = new double[numSets];
		double[] variance = new double[numSets];
		double[] logVariance = new double[numSets];
		long initTime; // For timings.

    StringBuffer str = new StringBuffer ("\n\n --------------------------");
		str.append (methodLabel + "\n  MC Variance : " + varMC + "\n\n");

		// Array-RQMC experiment with each pointSet.
		for (int i = 0; i < numSets; ++i) {
			initTime = System.currentTimeMillis();
			n = pointSets[i].getNumPoints();
			str.append ("n = " + n + "\n");
			simulReplicatesArrayRQMC (pointSets[i], rand, sort, sortCoordPts, 
          numSteps, m, statPerf);
			logn[i] = Num.log2(n);
			variance[i] = statPerf.variance();
			logVariance[i] = Num.log2(variance[i]);
			str.append ("  Average = " + statPerf.average() + "\n");
			str.append ("  VRF =  " + varMC / (n * variance[i]) + "\n");
			str.append (formatTime ((System.currentTimeMillis() - initTime) / 1000.) + "\n");
		}
		// Estimate regression slope and print plot and overall results.
	  double regSlope = slope (logn, logVariance, numSets);
		str.append ("Regression slope (log) for variance = " + regSlope  + "\n\n");

		// Print plot and overall results in files.
    if (filenamePlot != null) 
		try {
			Writer file = new FileWriter (filenamePlot + ".tex");
		  XYLineChart chart = new XYLineChart();
		  // ("title", "$log_2(n)$", "$log_2 Var[hat mu_{rqmc,s,n}]$");
		  chart.add (logn, logVariance);
			file.write (chart.toLatex(12, 8));
			file.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return str.toString();
   }

   /**
    * Sorts the chains that have not stopped yet using the stored
    * @ref umontreal.ssj.util.MultiDimSort. All the stopped chains are
    * placed at the end, then the chains that have not stopped are sorted.
    */
   public void sortNotStoppedChains (MultiDimSort sort) {
      int j = n - 1;
      int i = 0;
      T mc;
      while (j >= 0 && chains[j].hasStopped())
         --j;
      while (i < n && !chains[i].hasStopped())
         ++i;
      while (i < j) {
         while (!chains[i].hasStopped())
            ++i;
         while (chains[j].hasStopped())
            --j;
         mc = chains[i];
         chains[i] = chains[j];
         chains[j] = mc;
      }
      sort.sort(chains, 0, i);
   }

   /**
    * Sorts the chains using the stored
    * @ref umontreal.ssj.util.MultiDimSort.
    */
   public void sortChains() {
      savedSort.sort(chains, 0, n);
   }


  // Takes time in seconds and formats it.
	public String formatTime (double time) {
		int second, hour, min, centieme;
		hour = (int) (time / 3600.0);
		if (hour > 0) {
			time -= ((double) hour * 3600.0);
		}
		min = (int) (time / 60.0);
		if (min > 0) {
			time -= ((double) min * 60.0);
		}
		second = (int) time;
		centieme = (int) (100.0 * (time - (double) second) + 0.5);
		return String.valueOf(hour) + ":" + min + ":" + second + "." + centieme;
	}

	// Compute slope of linear regression of y on x, using only first n observations.
	public double slope(double[] x, double[] y, int n) {
		if (n < 2) {
			return 0.0;
		} else {
			double[] x2 = new double[n], y2 = new double[n];
			for (int i = 0; i < n; ++i) {
				x2[i] = x[i];
				y2[i] = y[i];
			}
			return LeastSquares.calcCoefficients (x2, y2, 1)[1];
		}
	}
}
