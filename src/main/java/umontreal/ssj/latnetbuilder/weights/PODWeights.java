package umontreal.ssj.latnetbuilder.weights;

import java.io.IOException;
import java.util.List;

import umontreal.ssj.mcqmctools.anova.CoordinateSet;
import umontreal.ssj.mcqmctools.anova.CoordinateSetLong;

/**
 * This class implements so-called product- and order dependent weights (POD weights). As the name suggests, POD weights are composed of @ref
 * OrderDependentWeights and @ref ProductWeights via @f$\gamma_{\mathfrak u} =
 * \Gamma_{|\mathfrak u|}\prod_{j\in\mathfrak u}\gamma_j@f$.
 * 
 * @author florian
 *
 */

public class PODWeights extends Weights {

	/**
	 * The product weight-component.
	 */
	ProductWeights productWeights;
	/**
	 * The order dependent weight-component.
	 */
	OrderDependentWeights orderDependentWeights;

	/**
	 * Constructs an instance of POD weights.
	 * @param pWeights The product weight-component.
	 * @param odWeights The order dependent weight-component.
	 */
	public PODWeights(ProductWeights pWeights, OrderDependentWeights odWeights) {
		this.productWeights = pWeights;
		this.orderDependentWeights = odWeights;
	}

	/**
	 * Constructor for given #productWeights and #orderDependentWeights is constructed by the default constructor for @ref OrderDependentWeights.
	 * @param pWeights The product weight-component.
	 */
	public PODWeights(ProductWeights pWeights) {
		this.productWeights = pWeights;
		this.orderDependentWeights = new OrderDependentWeights();
	}

	/**
	 * Constructor for given #orderDependentWeights and #productWeights  is constructed by the default constructor for @ref ProductWeights.
	 * @param odWeights The order dependent weight-component.
	 */
	public PODWeights(OrderDependentWeights odWeights) {
		this.productWeights = new ProductWeights();
		this.orderDependentWeights = odWeights;
	}
	
	
	/**
	 * Default constructor.
	 */
	public PODWeights() {
		this.productWeights = new ProductWeights();
		this.orderDependentWeights = new OrderDependentWeights();
	}

	/**
	 * Getter for the product weight-component.
	 * @return The product weight-component.
	 */
	public ProductWeights getProductWeights() {
		return productWeights;
	}
	
	/**
	 * Setter for the product weight-component.
	 */
	public void setProductWeights(ProductWeights productWeights) {
		this.productWeights = productWeights;
	}

	/**
	 * Getter for the order dependent weight-component.
	 * @return The order dependent weight-component.
	 */
	public OrderDependentWeights getOrderDependentWeights() {
		return orderDependentWeights;
	}

	/**
	 * Setter for the order dependent weight-component.
	 */
	public void setOrderDependentWeights(OrderDependentWeights orderDependentWeights) {
		this.orderDependentWeights = orderDependentWeights;
	}

	/**
	 * Adds the order dependent weight 'weight'.
	 * @param weight Order dependent weight to be added.
	 */
	public void addOrderDependentWeight(SingletonWeightComparable<Integer> weight) {
		orderDependentWeights.add(weight);
	}

	/**
	 * Adds an order dependent weight by passing the order and the value of the weight.
	 * @param ord The order of the weight to be added.
	 * @param weight The value of the weight to be added.
	 */
	public void addOrderDependentWeight(int ord, double weight) {
		orderDependentWeights.add(ord, weight);
	}
	
	/**
	 * Adds the product weight 'weight'.
	 * @param weight Product weight to be added.
	 */
	public void addProductWeight(SingletonWeightComparable<Integer> weight) {
		productWeights.add(weight);
	}
	
	/**
	 * Adds a product weight by passing the coordinate index and the value of the weight.
	 * @param index The coordinate index of the weight to be added.
	 * @param weight The value of the weight to be added.
	 */
	public void addProductWeight(int index, double weight) {
		productWeights.add(index, weight);
	}


	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer("");
		sb.append("POD weights:\n");
		sb.append(productWeights.toString() + "\n");
		sb.append(orderDependentWeights.toString());
		return sb.toString();
	}

	@Override
	public String toLatNetBuilder() {
		StringBuffer sb = new StringBuffer("");
		sb.append("POD:" + orderDependentWeights.getDefaultWeight());
		if (orderDependentWeights.weights.size() > 0) {
			sb.append(":");
			sb.append(orderDependentWeights.printBody());
		}
		sb.append(":" + productWeights.getDefaultWeight());
		if (productWeights.weights.size() > 0) {
			sb.append(":");
			sb.append(productWeights.printBody());
		}
		return sb.toString();
	}

}
