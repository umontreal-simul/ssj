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
