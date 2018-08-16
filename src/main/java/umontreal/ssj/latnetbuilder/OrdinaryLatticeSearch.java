/*
 * Class:        OrdinaryLatticeSearch
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

import umontreal.ssj.latnetbuilder.Search;

import umontreal.ssj.hups.Rank1Lattice;

import java.util.ArrayList;

/**
 * Class for the search of good rank-1 ordinary lattice rules using LatNet Builder.
 */
public class OrdinaryLatticeSearch extends Search{
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public Rank1Lattice search() throws RuntimeException
	{
		ArrayList<String> res = executeCommandLine();
		int numPoints = Integer.parseInt(res.get(1).split("  //")[0]); 
		int dimension = Integer.parseInt(res.get(2).split("  //")[0]);
		int[] genVec = new int[dimension];
		for(int coord = 0; coord < dimension; ++coord){
			genVec[coord] = Integer.parseInt(res.get(5 + coord).split("  //")[0]);
		}
		this.merit = Double.parseDouble(res.get(5 + dimension).split("  //")[0]);
		this.time = Double.parseDouble(res.get(6 + dimension).split("  //")[0]);
		this.successful = true;
		return new Rank1Lattice(numPoints, genVec, dimension);
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public String pointSetType()
	{
		return "lattice";
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public int interlacing()
	{
		return 1;
	}
	
	/**
	 * {@inheritDoc}
	 */
	@Override
	public String construction(){
		return "ordinary";
	}
}
