package umontreal.ssj.latnetbuilder.weights;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Implements comparable @ref Weights.
 * 
 * @author florian
 *
 * @param <T> indicates the data type of the index of the weight (e.g @ref
 *            Integer). Needs to be comparable.
 */


public abstract class WeightsComparable<T extends Comparable<T>> extends Weights<T> {

	protected ArrayList<SingletonWeightComparable<T>> weights; // actual weights
	boolean sorted = false;

	/**
	 * Constructs weights from a list of @ref SingletonWeightsComparable.
	 * 
	 * @param w list of weights with index type 'T'.
	 */
	public WeightsComparable(List<SingletonWeightComparable<T>> w) {
		weights = new ArrayList<SingletonWeightComparable<T>>(w);
	}

	/**
	 * Initializes an empty list of comparable weights.
	 */
	public WeightsComparable() {
		weights = new ArrayList<SingletonWeightComparable<T>>();
	}

	/**
	 * Returns the current comparable weights.
	 * 
	 * @return the list of comparable weights.
	 */
	// TODO: changing return type with getWeights() does not work. do we really need
	// this?
	public ArrayList<SingletonWeightComparable<T>> getComparableWeights() {
		return weights;
	}
	


	/**
	 * Adds a new comparable weight to the list. In case the weight for the
	 * respective index had already been set, it is overwritten.
	 * 
	 * @param singletonWeight comparable weight to be added.
	 */
	public void add(SingletonWeightComparable<T> singletonWeight) {
		boolean added = false;
		for (SingletonWeightComparable<T> w : weights) {
			if (w.getIndex() == singletonWeight.getIndex() && (!added)) {
				weights.set(weights.indexOf(w), singletonWeight);
				added = true;
			}
		}
		if (!added)
			weights.add(singletonWeight);
		
		sorted = false;
	}
	


	/**
	 * Sorts the weights w.r.t. the ordering defined on 'T', i.e., on the indices.
	 */
	public void sort() {
		Collections.sort(weights);
		sorted = true;
	}

	/**
	 * Adds a new weight with index 'index' and weight 'weight' or overwrites it, if
	 * the index already exists in the list.
	 * 
	 * @param index  index of the weight to be added.
	 * @param weight value of the weight to be added.
	 */
	public void add(T index, double weight) {
		add(new SingletonWeightComparable<T>(index, weight));
	}


}
