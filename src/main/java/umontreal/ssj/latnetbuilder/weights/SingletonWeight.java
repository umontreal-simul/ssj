package umontreal.ssj.latnetbuilder.weights;

/**
 * Implements a single instance of a weight for search-algorithms for quasi-Monte Carlo point sets
 *  as, e.g., those implemented in *LatNet Builder*. The template parameter 'T' indicates the data type of the 
 *  underlying index. Common choices are integers (e.g. the order for @ref OrderDependentWeights) or 
 *  @ref CoordinateSet (e.g. the projection for @ref ProjectionDependentWeights), etc.
 *  
 * @author florian
 *
 *@param <T> indicates the data type of the index of the weight (e.g. @ref CoordinateSet, @ref Integer,...).
 */

public class SingletonWeight<T> {
	
	
	protected T index;
	protected double weight;
	
	/**
	 * Constructs a weight with given 'index' and 'weight'.
	 * @param ind the desired index.
	 * @param w the desired value for the weight.
	 */
	SingletonWeight(T ind, double w){
		index = ind;
		weight = w;
	}
	
	/**
	 * Returns the current 'index'.
	 * @return the index of the weight.
	 */
	public T getIndex() {
		return index;
	}
	
	/**
	 * Sets the index of the weight to 'ind'.
	 * @param ind the desired value for the 'index'.
	 */
	public void setIndex(T ind) {
		index = ind;
	}
	
	/**
	 * Returns the current value of the weight.
	 * @return the value of the weight.
	 */
	public double getWeight() {
		return weight;
	}
	
	/**
	 * Sets the value of the weight to 'w'. 
	 * @param w desired value for the weight
	 */
	public void setWeight(double w) {
		weight = w;
	}
}
