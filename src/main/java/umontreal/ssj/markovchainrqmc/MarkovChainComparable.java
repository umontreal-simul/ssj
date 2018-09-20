package umontreal.ssj.markovchainrqmc;
import umontreal.ssj.util.sort.MultiDimComparable;

/**
 * A subclass of `MarkovChain` for which there is a total ordering between the
 * states, induced by the implementation of the
 * @ref umontreal.ssj.util.MultiDimComparable interface. A list or
 * array of `MarkovChainComparable` objects can be sorted according to their states at a
 * given step, using a  @ref umontreal.ssj.util.MultiDimComparator and an
 * external sorting method.
 *
 * The method `compareTo(MarkovChainComparable m, int i)` must return a
 * negative integer, zero, or a positive integer, and `this` Markov chain
 * is considered smaller than, equal to, or greater than the other chain `m` in dimension
 * `i`. The meaning of smaller, equal or larger is very flexible as well as
 * the meaning of dimension.
 *
 * Concrete subclasses must implement `compareTo` and must initialize
 * the variable `stateDim` in the constructor to the value of the state
 * dimension. It can be understood as the largest integer `j` for which the
 * method `compareTo(m,j)` can be called. They must also implement abstract
 * methods of `MarkovChain`.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class MarkovChainComparable extends MarkovChain
          implements MultiDimComparable<MarkovChainComparable>
{

   protected int stateDim;         // Dimension of the state

   /**
    * Returns the dimension of the state.
    */
   public int getStateDimension() {
     return stateDim;
  }

}
