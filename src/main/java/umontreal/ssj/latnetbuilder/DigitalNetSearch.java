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
import java.lang.Math;

/**
 * Class for the search of good digital nets using LatNet Builder.
 */
public class DigitalNetSearch extends Search {

	/**
	 * Class for the construction od digital nets.
	 */
	protected class DigitalNetBase2FromLatNetBuilder extends DigitalNetBase2 {
		
		public DigitalNetBase2FromLatNetBuilder (int numRows, int numCols, int dim, int[] matrices) {
			this.numCols = numCols;
			this.numRows = Math.min(numRows, MAXBITS);
			this.numPoints = 1 << this.numCols;
			this.dim = dim;
			this.genMat = matrices;
			this.outDigits = MAXBITS;
			// this.outDigits = this.numRows;
			this.normFactor = 1.0 / ( (double) (1L << this.outDigits) );
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

		dimension = dimension / interlacing;

		int[] genMat = new int[dimension*numCols];
		int trueNumRows = Math.min(31, numRows*interlacing);
		for(int coord = 0; coord < dimension; ++coord){
			for(int col = 0; col < numCols; ++col){
				genMat[coord * numCols + col] = 0;
				for(int row = 0; row < trueNumRows; ++row){
					// genMat[coord * numCols + col] += (1 << (trueNumRows - 1 - row)) * mats[coord*interlacing + row % interlacing][row/interlacing][col];
					genMat[coord * numCols + col] += (1 << (31 - 1 - row)) * mats[coord*interlacing + row % interlacing][row/interlacing][col];
				}
			}
		}
		this.merit = Double.parseDouble(res.get(res.size() - 2).split("  //")[0]);
		this.time = Double.parseDouble(res.get(res.size() - 1).split("  //")[0]);
		this.successful = true;
		return new DigitalNetBase2FromLatNetBuilder(trueNumRows, numCols, dimension, genMat);
	}
	
	/**
	 * Offset for the parsing of generating matrices.
	 */
	private int offsetForParsingGeneratingMatrix(int dimension){
		if (construction=="sobol"){
			return 7 + dimension;
		}
		else if (construction == "explicit"){
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
	public int interlacing()
	{
		return interlacing;
	}

	/**
	 * Sets the interlacing factor of the searched digital net.
	 * @param interlacing Interlacing factor.
	 */
	public void setInterlacing(int interlacing){
		this.interlacing = interlacing;
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
