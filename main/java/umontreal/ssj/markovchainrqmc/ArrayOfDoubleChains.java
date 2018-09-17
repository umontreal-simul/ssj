package umontreal.ssj.markovchainrqmc;

import umontreal.ssj.stat.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.util.*;
import umontreal.ssj.util.sort.*;
import java.util.Arrays;

/**
 * Similar to  @ref ArrayOfComparableChains, except that instead of working
 * with @f$n@f$ clones of a  @ref MarkovChain, we use a *single*
 * @ref MarkovChainDouble object for all the chains. The states of the chains
 * are maintained in an array of real numbers (<tt>double</tt>) and the
 * MarkovChainDouble.nextStepDouble method is used to advance each chain by
 * one step. The performance measure is assumed to be additive over all steps
 * of all copies of the chain. The sum is cumulated in a *single* accumulator
 * for all copies of the chain, updated at each step of each copy.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class ArrayOfDoubleChains extends ArrayOfComparableChains {

   protected MarkovChainDouble baseChain; // Single object for all chains.
   protected double[] state;              // States of the chains.
   protected double[] perfState;          // performance mesure of the states

   /**
    * Creates a virtual array for the chain <tt>baseChain</tt>. The method
    * <tt>makeCopies</tt> must be called to make the copies.
    */
   public ArrayOfDoubleChains (MarkovChainDouble baseChain,
                               PointSetRandomization rand) {
        super (baseChain, rand, new OneDimSort(0));
        this.baseChain = baseChain;
   }

   /**
    * Same as
    * {@link #ArrayOfDoubleChains(MarkovChainDouble,PointSetRandomization)
    * ArrayOfDoubleChains(baseChain, new RandomShift(new MRG32k3a()))}.
    */
   public ArrayOfDoubleChains (MarkovChainDouble baseChain) {
        this(baseChain, new RandomShift(new MRG32k3a()));
   }

   /**
    * Creates the vector of states for <tt>n</tt> copies of the base
    * chain.
    */
   public void makeCopies (int n) {
      state = new double[n];
      perfState = new double[n];
      this.n = n;
   }

   /**
    * Sets the states of the `n` copies of the base chain to `S`.
    */
   public void setStatesDouble (double[] S) {
      if(S.length!=n){
         this.n = S.length;
         System.out.println("WARNING : number of chains modified"+
         " to fit size of S in setStatesDouble(S)");
      }
      for (int i = 0; i < n; i++) {
          state[i] = S[i];
      }
   }

   /**
    * Returns the array containing the states of the `n` chains.
    */
   public double[] getStatesDouble() {
      return state;
   }

   /**
    * Initializes the states of the `n` copies of the base chain.
    */
   public void initStatesDouble() {
      Arrays.fill(state, baseChain.initialStateDouble());
      Arrays.fill(perfState, 0.0);
   }

   /**
    * Simulate one step for the `n` copies of the base chain, assuming
    * that we are at step `step`. The points are randomized before the
    * simulation using the stored
    * @ref umontreal.ssj.hups.PointSetRandomization. The dimension of `p`
    * must be at least as large as the number of uniforms required to
    * simulate one step of the chain. Returns true if and only if all
    * chains have stopped. Compute and stores the performances.
    */
   public boolean simulOneStepArrayRQMC (int step, PointSet p) {
      boolean allStopped = true;
      p.randomize(randomization);        // Randomize point set.
      PointSetIterator stream = p.iterator ();
      stream.resetStartStream (); // Beginning of pt set
      for (int i = 0; i < n; i++) {
          if(state[i] == Double.POSITIVE_INFINITY) continue;
          state[i] = baseChain.nextStepDouble (step, state[i], stream);
          stream.resetNextSubstream ();
          if(state[i] == Double.POSITIVE_INFINITY){
             perfState[i] = baseChain.getPerformance();
          }else{
             perfState[i] = baseChain.getPerformanceDouble(state[i],step);
          }
          allStopped = allStopped && state[i] == Double.POSITIVE_INFINITY;
      }
      return allStopped;
   }

   /**
    * Simulates the @f$n@f$ copies of the chain, <tt>numSteps</tt> steps
    * for each copy, using point set <tt>p</tt>, where @f$n@f$ is the
    * current number of copies of the chain and is *assumed* to equal the
    * number of points in <tt>p</tt>. At each step, the points are
    * randomized using the stored
    * @ref umontreal.ssj.hups.PointSetRandomization. The dimension of
    * <tt>p</tt> must be at least as large as the number of uniforms
    * required to simulate one step of the chain. Returns the average
    * performance per run.
    */
   public double simulArrayRQMC (PointSet p, int numSteps) {
      boolean allStopped = false;
      initStatesDouble();

      for (int step = 0; step < numSteps && !allStopped; step++) {
         sortChains();
         allStopped = simulOneStepArrayRQMC(step, p);

      }
      return calcMeanPerf();
   }

   /**
    * Computes and returns the mean performance of the @f$n@f$ chains.
    */
   public double calcMeanPerf() {
        double sumPerf = 0.0;                     // Sum of performances.
        for (int i=0; i<n; ++i) {
            sumPerf += perfState[i];
        }
        return sumPerf/n;
    }

   /**
    * Sorts the arrays containing the states of the @f$n@f$ chains.
    */
   public void sortChains() {
        Arrays.sort(state);
   }

   /**
    * Creates a String with the states.
    */
   public String toString() {
      StringBuffer sb = new StringBuffer(baseChain.toString());
      sb.append("***************************************************************"+PrintfFormat.NEWLINE);
      for(int i=0;i<n;++i)
         sb.append (" ; " + PrintfFormat.g(15, 6, state[i]));
      sb.append("PrintfFormat.NEWLINE");

      return sb.toString();
   }
} 
