package umontreal.ssj.util.sort.florian;

import java.util.Arrays;
import java.util.Comparator;

import umontreal.ssj.util.sort.MultiDimSort;

/**
 * This abstract class takes multi-dimensional objects such as arrays or
 * #MultiDim's, maps them to the real numbers via a #scoreFunction @f$h@f$ and
 * then sorts the objects according to their value under @f$h@f$.
 * 
 * The score function is individual for each setting and is thus implemented as
 * an abstract method.
 * 
 * @author florian
 *
 */
public abstract class MultDimToOneDimSort implements MultiDimSort<MultiDim> {

	protected int dimension;

	private double[][] indexForSort;

	private static class DoubleIndexComparator2 implements Comparator<double[]> {

		public int compare(double[] p1, double[] p2) {
			if (p1[1] > p2[1])
				return 1;
			else if (p1[1] < p2[1])
				return -1;
			else
				return 0;
		}
	}

	
	public void sort(MultiDim[] a, int iMin, int iMax) {
//		if (iMin == iMax)
//			return;
		double b[][] = new double[iMax][dimension];
		for (int i = iMin; i < iMax; ++i) {
			b[i] = a[i].getPoint();
//			for(int j = 0; j < a[i].getState().length; j++)
//			b[i][j] = a[i].getState()[j];
		}
		sort(b, iMin, iMax);

		// Now use indexForSort to sort a.
		// We do not want to clone all the objects in a,
		// but only the array of pointers.
		MultiDim[] aclone = a.clone(); // new Object[iMax];
//		MultiDim[] aclone = new MultiDim[a.length];
		for (int i = iMin; i < iMax; ++i)
			a[i] = aclone[(int) indexForSort[i][0]];

	}

	/**
	 * Function @f$h@f$ which maps the array \a b to the  real numbers.
	 * @param b
	 * @return
	 */
	public abstract double scoreFunction(double[] b);

	public void sort(MultiDim[] a) {
		sort(a, 0, a.length);
	}


	public int dimension() {
		return dimension;
	}

	public String toString() {
		return "MultDimToOneDimSort";
	}

	/**
	 * Sorts the `index` table by its second coordinate.
	 */
	public static void sortIndexOfDouble2(double[][] index, int iMin, int iMax) {
		Arrays.sort(index, iMin, iMax, new DoubleIndexComparator2());
	}

	public void sort(double[][] a, int iMin, int iMax) {
//		if (iMin + 1 == iMax)
//			return;
		indexForSort = new double[iMax][dimension];

		for (int i = iMin; i < iMax; ++i) {
			indexForSort[i][0] = i;
			indexForSort[i][1] = scoreFunction(a[i]);
		}
		Arrays.sort(indexForSort, iMin, iMax, new DoubleIndexComparator2());


		// Now use indexForSort to sort a.
		// We do not want to clone all the objects in a,
		// but only the array of pointers.
		double[][] aclone = a.clone(); 
		for (int i = iMin; i < iMax; ++i) {
			a[i] = aclone[(int) indexForSort[i][0]];
		}
	}

	public void sort(double[][] a) {
		sort(a, 0, a.length);

	}

}