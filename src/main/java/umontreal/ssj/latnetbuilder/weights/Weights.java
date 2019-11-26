package umontreal.ssj.latnetbuilder.weights;

import java.util.ArrayList;
import java.util.List;

/**
 * This abstract class implements generic weights that can be used to search for
 * good quasi-Monte Carlo point sets such as @ref Rank1Lattice, @ref DigitalNet,
 * etc. via search tools like *LatNet Builder.* The template parameter 'T'
 * indicates the data type of the underlying index. Common choices are integers
 * (e.g. the order for @ref OrderDependentWeights, @ref CoordinateSet for @ref
 * ProjectionDependentWeights, etc.). The abstract method #toLatNetBuilder is
 * used to communicate with the commandline tool of LatNet Builder.
 * 
 * @author florian
 *
 * @param <T> indicates the data type of the index of the weight (e.g. @ref
 *            CoordinateSet, @ref Integer,...).
 */
public abstract class Weights<T> {
	/**
	 * Weight to be used for indices that were not explicitly set. Its default value is zero.
	 */
	protected double defaultWeight = 0.0;
	protected ArrayList<SingletonWeight<T>> weights; // actual weights

	/**
	 * Constructs weights from a list of @ref SingletonWeights.
	 * 
	 * @param w list of weights with index type 'T'.
	 */
	public Weights(List<SingletonWeight<T>> w) {
		weights = new ArrayList<SingletonWeight<T>>(w);
	}

	/**
	 * Initializes an empty list of weights.
	 */
	public Weights() {
		weights = new ArrayList<SingletonWeight<T>>();
	}

	/**
	 * Returns the current weights as a list.
	 * 
	 * @return the list of weights.
	 */
	public ArrayList<SingletonWeight<T>> getWeights() {
		return weights;
	}

	/**
	 * Sets 'dWeight' as the current 'defaultWeight'.
	 * 
	 * @param dWeight desired default weight.
	 */
	public void setDefaultWeight(double dWeight) {
		defaultWeight = dWeight;
	}

	/**
	 * Returns the current 'defaultWeight'.
	 * 
	 * @return the value of the default weight.
	 */
	public double getDefaultWeight() {
		return defaultWeight;
	}

	/**
	 * Adds a new weight. In case the weight for the respective index had already
	 * been set, it is overwritten.
	 * 
	 * @param singletonWeight weight to be added.
	 */
	public void add(SingletonWeight<T> singletonWeight) {
		boolean added = false;
		for (SingletonWeight<T> w : weights) {
			if (w.getIndex() == singletonWeight.getIndex() && (!added)) {
				weights.set(weights.indexOf(w), singletonWeight);
				added = true;
			}
		}
		if (!added)
			weights.add(singletonWeight);
	}

	/**
	 * Adds a new weight with index 'index' and weight 'weight' or overwrites it, if
	 * the index already exists in the list.
	 * 
	 * @param index  index of the weight to be added.
	 * @param weight value of the weight to be added.
	 */
	public void add(T index, double weight) {
		add(new SingletonWeight<T>(index, weight));
	}

	/**
	 * Basic formatted string-output.
	 */
	public String toString() {
		StringBuffer sb = new StringBuffer("");
		sb.append("Weights [default = " + getDefaultWeight() + "]\n");

		sb.append("[");
		for (SingletonWeight<T> w : weights)
			sb.append(w.getWeight() + ",");
		sb.deleteCharAt(sb.length() - 1);
		sb.append("]\n");
		return sb.toString();
	}

	/*
	 * Methods to be implemented for non-abstract class
	 */
	/**
	 * Provides a String that can be interpreted by the command line interface of
	 * LatNetBuilder.
	 * 
	 * @return String for LatNetBuilder.
	 */
	public abstract String toLatNetBuilder();

}
