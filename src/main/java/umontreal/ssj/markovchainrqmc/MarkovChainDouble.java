package umontreal.ssj.markovchainrqmc;
import umontreal.ssj.rng.*;
import umontreal.ssj.stat.Tally;

/**
 * A special kind of Markov chain whose state space is a subset of the real
 * numbers. The state is stored in the double variable `state`.
 *
 * The "Double" version of the methods from Markov Chain now returns the
 * state. This makes it possible to simulate several copies of this chain in
 * parallel without cloning and without maintaining the state of the chain in
 * a local variable. The states can be maintained in an external array and at
 * each step, one passes the current state to the method
 * <tt>nextStepDouble</tt>, which returns the next state. This is exploited
 * in the implementation of <tt>ArrayOfDoubleChains</tt>.
 *
 * The methods  #initialState,  #nextStep, <tt>getPerformance</tt> and
 * #compareTo, which are abstract in  @ref MarkovChainComparable, all have a
 * default implementation here, so the methods of Markov Chain still work.
 *
 * On the other hand, the abstract methods specified in the present class do
 * not necessarily update local variables.
 *
 * Abstract method `compareTo` from class `MarkovChainComparable` is
 * implemented by comparing variable `state` of the chains. And method
 * `dimension` now returns 1, there’s no need to set variable `stateDim`.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class MarkovChainDouble extends MarkovChainComparable {

   protected double state = 0;             // State of this chain.
   protected int step = 0;                 // Current step.
   protected double perf = 0;

   /**
    * @name Abstract Methods
    * @{
    */

   /**
    * Returns the initial (deterministic) state.
    */
   public abstract double initialStateDouble();

   public void initialState() {
        step = 0;
        state = initialStateDouble();
        stopped = false;
   }

   /**
    * Simulates one step of the chain, from state `s`, using `stream` for
    * the randomness, assuming we are at step `step`. If the chain stops,
    * this method must compute the performance and save it to variable
    * `perf`, then return <tt>Double.POSITIVE_INFINITY</tt>. If not, it
    * returns the new state.
    */
   public abstract double nextStepDouble (int step, double s,
                                          RandomStream stream);

   public void nextStep (RandomStream stream) {
        state = nextStepDouble (step, state, stream);
        step++;
   }

   /**
    * Returns the performance measure associated with state `state`, which
    * may depend on the number of steps `numsteps`.
    */
   public abstract double getPerformanceDouble (double state, int numSteps) ;

   /**
    * @}
    */

   /**
    * @name Other Methods
    * @{
    */

   /**
    * Returns the performance mesure associated with current state, which
    * may depend on the number of steps `numsteps`.
    */
   public double getPerformance (int numSteps) {
      return getPerformanceDouble(state,numSteps);
   }

   /**
    * Returns the value of `perf` which is computed when a chain stops.
    */
   public double getPerformance() {
       return perf;
    }

   /**
    * Indicates if the chain has stopped.
    */
   public boolean hasStopped() {
        return state == Double.POSITIVE_INFINITY;
    }


   public int getStateDimension(){
      return 1;
   }

   public int compareTo (MarkovChainComparable other, int i) {
       double os = ((MarkovChainDouble)other).state;
       return (state < os ? -1 : (state > os ? 1 : 0));
   }

/**
 * After invoking  #initialStateDouble, starts a new simulation run,
 * simulates <tt>numSteps</tt> steps of the Markov chain using the given
 * <tt>stream</tt>, and returns the final state. The  #simulSteps method in
 * `MarkovChain` does the same, but returns nothing.
 */
   public double simulStepsDouble (int numSteps, RandomStream stream) {
      initialState();
      for(step = 0; step<numSteps && !hasStopped(); ++step){
         state = nextStepDouble(step,state,stream);
      }
      return state;
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
            simulStepsDouble (numSteps, stream);
            statRuns.add (getPerformanceDouble(state,numSteps));
            stream.resetNextSubstream ();
        }
    }

 }

/**
 * @}
 */