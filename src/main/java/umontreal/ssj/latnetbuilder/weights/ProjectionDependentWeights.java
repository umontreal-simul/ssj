package umontreal.ssj.latnetbuilder.weights;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import umontreal.ssj.mcqmctools.anova.CoordinateSet;
import umontreal.ssj.mcqmctools.anova.CoordinateSetLong;

/**
 * Class implementing *projection dependent weights.* In this case, weights are  assigned to each projection individually. With a view to using this
 * class together with *LatticeBuilder* it is also possible to set weights for the entire collections of projections of the same order via 
 * @ref OrderDependentWeights.
 * @author florian
 *
 */

public class ProjectionDependentWeights extends Weights<CoordinateSet> {

	public OrderDependentWeights orderDependentWeights; 
	protected String fileDir = "";
	protected String fileName = "";

	/**
	 * Sets the directory to which an output-file can be generated.
	 * @param dir path to an output directory.
	 */
	public void setFileDir(String dir) {
		 fileDir = dir;
	}
	
	/**
	 * Returns the directory to which an output-file can be generated.
	 * @return path to an output directory.
	 */
	public String getFileDir() {
		return fileDir;
	}
	
	/**
	 * Sets the name of the file to which an output can be generated.
	 * @param name name for output-file.
	 */
	public void setFileName(String name) {
		fileName = name;
	}
	
	/**
	 * Returns the name of the file to which an output can be generated.
	 * @return name for output-file.
	 */
	public String getFileName() {
		return fileName;
	}
	
	/**
	 * Instantiates the parameter #orderDependentWeights.
	 */
	public void addOrderDependentWeight() {
		orderDependentWeights = new OrderDependentWeights();
	}
	
	/**
	 * Adds an order dependent weight to #orderDependentWeights. Duplicate assignments to the same order are taken care of.
	 * @param weight the comparable 'SingletonWeight' that is to be added.
	 */
	public void addOrderDependentWeight(SingletonWeightComparable<Integer> weight) {
		orderDependentWeights.add(weight);
	}
	
	/**
	 * Adds an order dependent weight to #orderDependentWeights. Duplicate assignments to the same order are taken care of.
	 * @param order the order of the projections for which the weight is assigned.
	 * @param weight the value of the weight that is to be added.
	 */
	public void addOrderDependentWeight(int order, double weight) {
		orderDependentWeights.add(order,weight);
	}
	
	/**
	 * Assigns an entire list of comparable 'SingleTonWeights' as order dependent weights to #orderDependentWeights.
	 * @param ordWeights list of singleton weights representing the order dependent weights.
	 */
	public void setOrderDependentWeights(ArrayList<SingletonWeightComparable<Integer>> ordWeights) {
		orderDependentWeights = new OrderDependentWeights(ordWeights);
	}

	/**
	 * Writes a file with name #fileName to the directory #fileDir containing information on the weights. This file is formatted to be further
	 * processed by *LatticeBuilder.** 
	 * @throws IOException
	 */
	public void write() throws IOException {
		FileWriter file = new FileWriter(fileDir + fileName);
		StringBuffer sb = new StringBuffer("");
	
		for(SingletonWeight<CoordinateSet> w : weights)
			sb.append(w.getIndex().toStringNoBraces() + ":\t" + w.getWeight() + "\n");
		
		for(SingletonWeightComparable<Integer> w : orderDependentWeights.getComparableWeights())
			sb.append("order " + w.getIndex() + ":\t" + w.getWeight() + "\n");
		sb.append("default:\t" + defaultWeight);
		file.write(sb.toString());
		file.close();
	}
	
	/**
	 * Creates a formatted output of the projection dependent weights.
	 * @return a formatted output of the projection dependent weights.
	 */
	@Override
	public String toString(){
		StringBuffer sb = new StringBuffer("");
		sb.append("Order dependent weights [default = " + getDefaultWeight() + "]\n");
		
		for(SingletonWeight<CoordinateSet> w : weights)
			sb.append(String.format("%20s: %1.8g%n", w.getIndex().toString(),w.getWeight()));

		for(SingletonWeightComparable<Integer> w : orderDependentWeights.getComparableWeights())
			sb.append(String.format("%20s: %1.8g%n","order " + w.getIndex(),w.getWeight()));

		return sb.toString();
	}
	
	public  String toLatNetBuilder(){
		String str = "file:" + fileDir + fileName;
		return str;
	};
	
	
	
}
