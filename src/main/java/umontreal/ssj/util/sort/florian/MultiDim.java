
package umontreal.ssj.util.sort.florian;

/**
 * This interface represents a point or array of @f$d@f$ dimensions
 * in @f$\mathbb{R}^d@f$. The value of the @f$j@f$th dimension can be accessed
 * with the method {@link #getCoordinate() getCoordinate(j)}.
 * 
 * @remark **Florian:** This is what I used for my experiments with chemical
 *         reaction networks. It is basically just another copy of @ref
 *         MultiDim01, but without the restriction to the unit cube. (I have not
 *         merged them to not mess with the current structure). It can probably
 *         be embedded more elegantly into the existing structure with other interfaces of this package, but
 *         there is not enough time for that now. 
 */
public interface MultiDim {

	/**
	 * This method returns the number dimensions of this point.
	 */
	public int dimension();

	/**
	 * Returns the @f$d@f$ coordinates of this point.
	 */
	public double[] getPoint();

	/**
	 * Returns the value of @f$j@f$th coordinate (or dimension).
	 */
	public double getCoordinate(int j);



}
