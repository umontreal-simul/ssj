package umontreal.ssj.stat.density;

public interface DensityEstimationModel {
	
	
	public double estimateIV(DensityEstimator de);
	public double estimateISB(DensityEstimator de);
	public double estimateMISE(DensityEstimator de);
}
