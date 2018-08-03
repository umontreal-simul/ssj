package umontreal.ssj.latnetbuilder;

import umontreal.ssj.latnetbuilder.Search;

import umontreal.ssj.hups.Rank1Lattice;

import java.util.ArrayList;

public class OrdinaryLatticeSearch extends Search{
	
	@Override
	public Rank1Lattice search()
	{
		try {
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
		} catch (RuntimeException e)
		{
			return null;
		}
	}
	
	@Override
	public String pointSetType()
	{
		return "lattice";
	}
	
	@Override
	public int interlacing()
	{
		return 1;
	}
	
	@Override
	public String construction(){
		return "ordinary";
	}
}
