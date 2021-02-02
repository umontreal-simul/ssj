/*
 * Class:        DigitalNetSearch
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

import umontreal.ssj.hups.DigitalNetBase2;
import java.util.ArrayList;

/**
 * Class for the search of good digital nets using LatNet Builder.
 */
public class DigitalNetSearch extends Search {

	/**
	 * Class for the construction od digital nets.
	 */
	protected class DigitalNetBase2FromLatNetBuilder extends DigitalNetBase2 {
		
		public DigitalNetBase2FromLatNetBuilder (int numRows, int numCols, int dim, int interlacing, int[][][] matrices) {
			if (numRows > MAXBITS || numCols > MAXBITS)
				throw new RuntimeException (String.format("SSJ cannot handle matrices with more than %d rows or columns.", MAXBITS));
			this.numCols = numCols;
			this.numRows = numRows;
			this.numPoints = 1 << this.numCols;
			this.dim = dim;
			this.outDigits = MAXBITS;
			this.normFactor = 1.0 / ( (double) (1L << this.outDigits) );
			this.interlacing = interlacing;

			generatorMatricesFromStandardFormat(matrices);
		}
	}
	
	String construction;
	int interlacing;
	
	/**
	 * Constructor.
	 * @param construction Type of construction (eg. sobol, explicit, polynomial, ...). 
	 */
	public DigitalNetSearch(String construction)
	{
		super();
		this.construction = construction;
		this.interlacing = 1;
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public DigitalNetBase2 search() throws RuntimeException
	{
		ArrayList<String> res = executeCommandLine();
		int numCols = Integer.parseInt(res.get(0).split("  //")[0]); 
		int numRows = Integer.parseInt(res.get(1).split("  //")[0]);
		int interlacing = Integer.parseInt(res.get(4).split("  //")[0]);
		int dimension = Integer.parseInt(res.get(3).split("  //")[0]);

		int[][][] mats = new int[dimension][numRows][numCols];
		for(int coord = 0; coord < dimension; ++coord){
			for(int row = 0; row < numRows; ++row){
				String[] tmp = res.get(coord*(numRows+1)+row+offsetForParsingGeneratingMatrix(dimension)).split(" ");
				for(int col = 0; col < numCols; ++col){
					mats[coord][row][col] = Integer.parseInt(tmp[col]);
				}
			}
		}

		this.merit = Double.parseDouble(res.get(res.size() - 2).split("  //")[0]);
		this.time = Double.parseDouble(res.get(res.size() - 1).split("  //")[0]);
		this.successful = true;

		return new DigitalNetBase2FromLatNetBuilder(numRows, numCols, dimension, interlacing, mats);
	}
	
	/**
	 * Offset for the parsing of generating matrices.
	 */
	private int offsetForParsingGeneratingMatrix(int dimension){
		if (construction.equals("sobol")){
			return 7 + dimension;
		}
		else if (construction.equals("explicit")){
			return 7;
		}
		else
		{
			return 8 + dimension;
		}
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public String pointSetType()
	{
		return "net";
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public String interlacing()
	{
		return String.valueOf(interlacing);
	}

	/**
	 * Sets the interlacing factor of the searched digital net.
	 * @param interlacing Interlacing factor.
	 */
	public void setInterlacing(String interlacing){
		this.interlacing = Integer.parseInt(interlacing);
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public String construction(){
		return construction;
	}
	
	/**
	 * Sets the construciton method of the searched digital net.
	 * @param construction Type of construction (eg. sobol, explicit, polynomial, ...).
	 */
	public void setConstruction(String construction){
		this.construction = construction;
	}
	
}
