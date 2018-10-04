/*
 * Class:        Search
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2018  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author Maxime Godin and Pierre Marion
 * @since August 2018
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.latnetbuilder;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.File;
import java.io.FileReader;

import java.util.List;
import java.util.ArrayList;
import umontreal.ssj.hups.PointSet;

/**
 * Abstract class for the search of highly uniform point sets with LatNet Builder.
 * See the <a href="http://umontreal-simul.github.io/latnetbuilder/de/db5/cmdtut_summary.html">Summary of Command-Line Options</a>
 * of LatNet Builder for details about the various fields of this class.
 */
abstract public class Search {

	String pathToOutputFolder;
	
	static String PATH_TO_LATNETBUILDER = "latnetbuilder"; // path to the latnetbuilder executable
	
	String sizeParameter;
	int dimension;
	boolean multilevel;
	String combiner;
	String explorationMethod;
	String figure;
	String normType;
	ArrayList<String> weights;
	ArrayList<String> filters;

	boolean successful;
	double merit;
	double time;
	
	/**
	 * Constructor.
	 */
	protected Search()
	{
		this.pathToOutputFolder = "latnetbuilder_results";
		this.multilevel = false;
		this.weights = new ArrayList<String>();
		this.filters = new ArrayList<String>();
		this.successful = false;
		this.merit = 0;
		this.time = 0;
	}
	
	@Override
	/**
	 * Formats the search parameters for printing.
	 */
	public String toString() {
		return "Point Set Type: " + pointSetType() + "\n" +
				"Construction method: " + construction() + "\n" +
				"Size parameter: " + sizeParameter + "\n" +
				"Multilevel: " + multilevel + "\n" +
				"Dimension: " + dimension + "\n" +
				"Interlacing: " + interlacing() + "\n" +
				"Exploration method: " + explorationMethod + "\n" +
				"Figure of merit: " + figure + "\n" +
				"Norm-type: " + normType + "\n" + 
				"Combiner: " + combiner + "\n" + 
				"Filters: " + filters + "\n" + 
				"Weights: " + weights + "\n" +
				"Output folder: " + pathToOutputFolder;
				
	}
	
	/**
	 * Returns the type of the searched point set.
	 */
	abstract public String pointSetType();
	
	/**
	 * Returns the interlacing factor of the search.
	 */
	abstract public int interlacing();
	
	
	abstract public String construction();
	
	/**
	 * Constructs the command-line for the LatNet Builder executable.
	 */
	private String constructCommandLine() {
		StringBuffer sb = new StringBuffer(PATH_TO_LATNETBUILDER + " -v 0");
		sb.append(" -t " +  pointSetType());
		sb.append(" -c " +  construction());
		sb.append(" -M " +  multilevel);
		sb.append(" -s " +  sizeParameter);
		sb.append(" -f " +  figure);
		sb.append(" -q " +  normType);
		sb.append(" -e " +  explorationMethod);
		sb.append(" -d " +  dimension);
		sb.append(" -i " +  interlacing());
		sb.append(" --output-folder " +  pathToOutputFolder);
		sb.append(" -w");
		for(String w : weights)
		{
			sb.append(" " + w);
		}
		if (filters.size() > 0)
			sb.append(" -F ");
		for(String f : filters)
		{
			sb.append(" " + f);
		}	
		if (multilevel)
		{
			sb.append(" -C " + combiner);
		}
		return sb.toString();
	}
	
	/**
	 * Executes the command-line and reads the content of the outputMachine.txt file.
	 */
	protected ArrayList<String> executeCommandLine() {
		String cmd = constructCommandLine();
		BufferedReader br = null;
		try {
			Process p = Runtime.getRuntime().exec(cmd);
			BufferedReader stderr = new BufferedReader(new InputStreamReader(p.getErrorStream()));
	
			try {
				p.waitFor();
				String res = "";
				if (p.exitValue() != 0) {
					while (stderr.ready())
						res = res + "\n" + stderr.readLine();
					throw new RuntimeException("LatNet Builder exited with status " + p.exitValue() + "\nCOMMAND LINE: " + cmd + res);
				}
			} catch (InterruptedException e) {
				throw new RuntimeException("LatNet Builder interrupted");
			}
			File f = new File(pathToOutputFolder, "outputMachine.txt");
			br = new BufferedReader (new FileReader(f));
			
			ArrayList<String> res = new ArrayList<String>();
			String line;
			while((line = br.readLine()) != null){
				res.add(line);
			}
			return res;
		} catch (IOException e) { 
			throw new RuntimeException("Error in the communication with LatNet Builder"); 
		}
		finally { 
			try {
				if (br != null){
					br.close();
				}
			} catch (IOException e) {
				throw new RuntimeException("Error in the communication with LatNet Builder");
			}
		}
	}
	
	/**
	 * Executes the search and returns the corresponding point set.
	 */
	abstract public PointSet search() throws RuntimeException;

	/**
	 * Returns the merit value of the point set.
	 */
	public double merit(){
		return this.merit;
	}

	/**
	 * Returns the elapsed CPU time taken for the search.
	 */
	public double time(){
		return this.time;
	}

	/**
	 * Returns a boolean indicating if the search was successful.
	 */
	public boolean successful(){
		return this.successful;
	}
	
	/**
	 * Sets the path to the latnetbuilder executable.
	 * @param path Path to the latnetbuilder executable.
	 */
	public void setPathToLatNetBuilder(String path) {
		PATH_TO_LATNETBUILDER = path;
	}
	
	/**
	 * Sets the dimension of the searched point set.
	 * @param dimension Dimension of the point set.
	 */
	public void setDimension(int dimension) {
		this.dimension = dimension;
	}
	
	/**
	 * Sets the size parameter of the point set
	 * @param sizeParameter Modulus for lattices and number of points for digital nets.
	 */
	public void setSizeParameter(String sizeParameter) {
		this.sizeParameter = sizeParameter;
	}
	
	/**
	 * Sets the search to the multilevel mode.
	 * @param multilevel Flag for the multilevel mode.
	 */
	public void setMultilevel(boolean multilevel) {
		this.multilevel = multilevel;
	}
	
	/**
	 * Sets the combiner for the merit in the multilevel mode case.
	 * @param combiner Combiner.
	 */
	public void setCombiner(String combiner) {
		this.combiner = combiner;
	}
	
	/**
	 * Sets the exploration method of the search.
	 * @param explorationMethod Exploration method.
	 */
	public void setExplorationMethod(String explorationMethod) {
		this.explorationMethod = explorationMethod;
	}
	
	/**
	 * Sets the figure of merit of the search.
	 * @param figure Figure of merit.
	 */
	public void setFigureOfMerit(String figure) {
		this.figure = figure;
	}
	
	/**
	 * Sets the norm-type for the figure of merit of the search.
	 * @param normType Norm-type.
	 */
	public void setNormType(String normType) {
		this.normType = normType;
	}
	
	/**
	 * Sets the weights to the figure of merit of the search.
	 * @param weights List of weights.
	 */
	public void setWeights(List<String> weights) {
		this.weights = new ArrayList<String>(weights);
	}
	
	/**
	 * Add a weight to the figure of merit of the search.
	 * @param weight Weight.
	 */
	public void addWeight(String weight) {
		this.weights.add(weight);
	}
	
	/**
	 * Sets the filters of the search.
	 * @param filters List of filters.
	 */
	public void setFilters(List<String> filters) {
		this.filters = new ArrayList<String>(filters);
	}

	/**
	 * Sets the path to the output folder.
	 */
	public void setPathToOutputFolder(String path) {
		this.pathToOutputFolder = path;
	}
}