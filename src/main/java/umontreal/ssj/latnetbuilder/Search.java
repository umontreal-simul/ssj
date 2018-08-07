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
 */
abstract public class Search {

	String pathToOutputFolder;
	
	static String PATH_TO_LATNETBUILDER = ""; // path to the directory containing the latnetbuilder executable
	
	String sizeParameter;
	int dimension; // dimension of the point set
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
	
	abstract public String pointSetType();
	
	abstract public int interlacing();
	
	abstract public String construction();
	
	private String constructCommandLine() {
		StringBuffer sb = new StringBuffer(new File(PATH_TO_LATNETBUILDER,"latnetbuilder").toString() +" -v 0");
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
	
	abstract public PointSet search() throws RuntimeException;

	public double merit(){
		return this.merit;
	}

	public double time(){
		return this.time;
	}

	public boolean successful(){
		return this.successful;
	}
	
	public void setPathToLatNetBuilder(String path) {
		PATH_TO_LATNETBUILDER = path;
	}
	
	public void setDimension(int dimension) {
		this.dimension = dimension;
	}
	
	public void setSizeParameter(String sizeParameter) {
		this.sizeParameter = sizeParameter;
	}
	
	public void setMultilevel(boolean multilevel) {
		this.multilevel = multilevel;
	}
	
	public void setCombiner(String combiner) {
		this.combiner = combiner;
	}
	
	public void setExplorationMethod(String explorationMethod) {
		this.explorationMethod = explorationMethod;
	}
	
	public void setFigureOfMerit(String figure) {
		this.figure = figure;
	}
	
	public void setNormType(String normType) {
		this.normType = normType;
	}
	
	public void setWeights(List<String> weights) {
		this.weights = new ArrayList<String>(weights);
	}
	
	public void addWeight(String weight) {
		this.weights.add(weight);
	}
	
	public void setFilters(List<String> filters) {
		this.filters = new ArrayList<String>(filters);
	}

	public void setPathToOutputFolder(String path) {
		this.pathToOutputFolder = path;
	}
}