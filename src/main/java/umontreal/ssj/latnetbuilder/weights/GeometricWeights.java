package umontreal.ssj.latnetbuilder.weights;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Class implementing a special case of @ref OrderDependentWeights,@f$\Gamma_{k}
 * = c \rho^{k-1} @f$ , where each weight is given by a #baseWeight @f$\rho@f$
 * raised to @f$k-1@f$, i.e., the power of the order of the projection it is
 * assigned to - 1. By default, @f$\rho=1@f$. It is also possible to set a
 * #truncationLevel, so that for all orders @f$k@f$ larger than this number each
 * weight equals the #defaultWeight.
 * 
 * @author florian
 *
 */

public class GeometricWeights extends OrderDependentWeights {

	protected double baseWeight = 1.0;
	protected int truncationLevel;
	protected double c;

	/**
	 * Constructs geometric weights with given #baseWeight 'base' and
	 * #truncationLevel 'truncation'.
	 * 
	 * @param base       desired value for @f$\rho@f$.
	 * @param truncation order up to which geometric weights are set. Higher-order
	 *                   projections will be assigned the #defaultWeight.
	 */

	public GeometricWeights(double base, int truncation, double c) {
		super();
		baseWeight = base;
		truncationLevel = truncation;
		this.c = c;
		setWeights();
	}

	/**
	 * Constructor that sets @f$c=1@f$.
	 * 
	 * @param base       desired value for @f$\rho@f$.
	 * @param truncation rder up to which geometric weights are set. Higher-order
	 *                   projections will be assigned the #defaultWeight.
	 */

	public GeometricWeights(double base, int truncation) {
		this(base, truncation, 1.0);
	}

	/**
	 * Default constructor.
	 */
	public GeometricWeights() {
		super();
		truncationLevel = 0;
	}

	/**
	 * Getter for the constant factor @f$c@f$.
	 * 
	 * @return the constant factor@f$c@f$
	 */
	public double getC() {
		return c;
	}

	/**
	 * Sets the constant factor @f$c@f$.
	 * 
	 * @param c the constant factor @f$c@f$
	 */
	public void setC(double c) {
		this.c = c;
	}

	/**
	 * Returns the #baseWeight.
	 * 
	 * @return the base value @f$\rho@f$.
	 */
	public double getBaseWeight() {
		return baseWeight;
	}

	/**
	 * Sets the #baseWeight @f$\rho@f$.
	 * 
	 * @param base desired value for @f$\rho@f$.
	 */
	public void setBaseWeight(double base) {
		baseWeight = base;
	}

	/**
	 * Returns the current #truncationLevel. Higher-order projections will be
	 * assigned the #defaultWeight.
	 * 
	 * @return truncation order up to which the weights are computed.
	 */
	public int getTruncationLevel() {
		return truncationLevel;
	}

	/**
	 * Sets the #truncationLevel. Higher-order projections will be assigned the
	 * #defaultWeight.
	 * 
	 * @param trLevel desired order up to which the weights are computed.
	 */
	public void setTruncationLevel(int trLevel) {
		truncationLevel = trLevel;
	}


	/**
	 * Computes the weights up to the order #truncationLevel and assigns them. 
	 */
	public void setWeights() {
		double w = 1.0;
		weights = new ArrayList<SingletonWeightComparable<Integer>>(truncationLevel);
		for (int order = 1; order <= truncationLevel; order++) {
			weights.add(order - 1, new SingletonWeightComparable<Integer>(order, c * w));
			w *= baseWeight;
		}
	}

	/**
	 * Creates a formatted output of the geometric order dependent weights ordered
	 * w.r.t. to their order @f$k@f$.
	 * 
	 * @return a formatted output of the geometric weights.
	 */
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer("");
		sb.append("Geometric order dependent weights [default = " + getDefaultWeight() + ", base = " + baseWeight
				+ "]:\n");
		if (weights.size() > 0)
			sb.append("[");
		sb.append(printBody());
		return sb.toString() + (weights.size() > 0 ? "]" : "");
	}



}
