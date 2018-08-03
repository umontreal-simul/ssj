package umontreal.ssj.latnetbuilder;

import umontreal.ssj.hups.DigitalNetBase2;

import java.util.ArrayList;

public class DigitalNetSearch extends Search {

	protected class DigitalNetBase2FromLatNetBuilder extends DigitalNetBase2{
		
		DigitalNetBase2FromLatNetBuilder(int numRows, int numCols, int dim, int[] matrices) {
			this.numCols = numCols;
			this.numRows = numRows;
			this.numPoints = 1 << this.numCols;
			this.dim = dim;
			this.genMat = matrices;
			this.outDigits = MAXBITS;
			this.normFactor = 1.0 / ( (double) (1L << this.numRows) );
		}
	}
	
	String construction;
	
	public DigitalNetSearch(String construction)
	{
		super();
		this.construction = construction;
	}
	
	@Override
	public DigitalNetBase2 search()
	{
		try {
			ArrayList<String> res = executeCommandLine();
			int numCols = Integer.parseInt(res.get(0).split("  //")[0]); 
			int numRows = Integer.parseInt(res.get(1).split("  //")[0]);
			int dimension = Integer.parseInt(res.get(3).split("  //")[0]);
			int[] genMat = new int[dimension*numCols];
			for(int coord = 0; coord < dimension; ++coord){
				int[][] mat = new int[numRows][numCols];
				for(int row = 0; row < numRows; ++row){
					String[] tmp = res.get(coord*(numRows+1)+row+offsetForParsingGeneratingMatrix(dimension)).split(" ");
					for(int col = 0; col < numCols; ++col){
						mat[row][col] = Integer.parseInt(tmp[col]);
					}
				}
				for(int col = 0; col < numCols; ++col){
					genMat[coord * numCols + col] = 0;
					for(int row = 0; row < numRows; ++row){
						genMat[coord * numCols + col] += (1 << row) * mat[numRows - 1- row][col];
					}
				}
			}
			this.merit = Double.parseDouble(res.get(res.size() - 2).split("  //")[0]);
			this.time = Double.parseDouble(res.get(res.size() - 1).split("  //")[0]);
			this.successful = true;
			return new DigitalNetBase2FromLatNetBuilder(numRows, numCols, dimension, genMat);
		} catch (RuntimeException e) {
			return null;
		}
	}
	
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
	
	@Override
	public String pointSetType()
	{
		return "net";
	}
	
	@Override
	public int interlacing()
	{
		return 1;
	}
	
	@Override
	public String construction(){
		return construction;
	}
	
	public void setConstruction(String construction){
		this.construction = construction;
	}
	
}
