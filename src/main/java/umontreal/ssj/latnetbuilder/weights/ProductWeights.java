package umontreal.ssj.latnetbuilder.weights;

import java.util.ArrayList;

/**
 * Class implementing *product weights.* Each weight is assigned to one
 * coordinate, i.e., @f$\gamma_{\mathfrak u} = \prod_{j\in\mathfrak u}\gamma_j@f$. Thus, the index of each
 * 'SingletonWeightComparable' is given by an integer (starting at 0 for the
 * first coordinate) and corresponds to one specific coordinate.
 * 
 * @author puchhamf
 *
 */

public class ProductWeights extends WeightsComparable<Integer> {


	/**
	 * Constructor for product weights from a list of @ref SingletonWeightComparable indexed by integers indicating the coordinate.
	 * @param weightList list of comparable weights to be set as the product weights.
	 */
	public ProductWeights(ArrayList<SingletonWeightComparable<Integer>> weightList) {
		super(weightList);
	}

	/**
	 * Default constructor.
	 */
	public ProductWeights() {
		super();
	}


	/**
	 * Sorts the weights and creates a rudimentary string containing the product
	 * weights separated by commas. Missing weights (i.e. a weight for a coordinates
	 * which has not been set and lies between two coordinates, whose weights are
	 * specified) are filled with 'defaultWeight'.
	 * 
	 * @return a string containing the values of the weights separated by commas.
	 */
	String printBody() {
		if (!sorted)
			sort();
		StringBuffer sb = new StringBuffer("");
		if (weights.size() > 0) {
			int index = 0;
			for (SingletonWeightComparable<Integer> w : weights) {
				while (index < w.getIndex()) {
					sb.append(getDefaultWeight() + ",");
					index++;
				}
				sb.append(w.getWeight() + ",");
				index++;
			}
			sb.deleteCharAt(sb.length() - 1);
		}
		return sb.toString();
	}

	/**
	 * Creates a formatted output of the product weights ordered w.r.t. to the
	 * coordinate they are assigned to. Missing weights (i.e. a weight for a coordinates
	 * which has not been set and lies between two coordinates, whose weights are
	 * specified) are filled with 'defaultWeight'.
	 * 
	 * @return a formatted output of the product weights.
	 */
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer("");
		sb.append("Product weights [default = " + getDefaultWeight() + "]:\n");
		if (weights.size() > 0)
			sb.append("[");
		sb.append(printBody());
		return sb.toString() + (weights.size() > 0 ? "]" : "");
	}

	/**
	 * Creates a string formatted for passing it to *LatNet Builder*.
	 * 
	 * @return a formatted string that can be processed by LatNet Builder.
	 */
	public String toLatNetBuilder() {
		StringBuffer sb = new StringBuffer("");
		sb.append("product:" + getDefaultWeight());
		if (weights.size() > 0) {
			sb.append(":");
			sb.append(printBody());
		}
//		sb.append(" -o " + weightPower + " ");

		return sb.toString();
	}

	
}