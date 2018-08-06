/*
 * Class:        PolynomialLatticeSearch
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

import umontreal.ssj.latnetbuilder.DigitalNetSearch;

public class PolynomialLatticeSearch extends DigitalNetSearch{

	String pointSetType;
	
	public PolynomialLatticeSearch(String pointSetType){
		super("polynomial");
		this.pointSetType = pointSetType;
	}
	
	@Override
	public String pointSetType(){
		return pointSetType;
	}
	
	public void changePointSetTypeView(String pointSetType)
	{
		this.pointSetType = pointSetType;
	}
	
	@Override
	public void setConstruction(String construction) {
		return;
	}
}
