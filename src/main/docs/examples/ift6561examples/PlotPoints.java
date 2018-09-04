package ift6561examples;
import umontreal.ssj.util.PrintfFormat;
import umontreal.ssj.hups.*;
import umontreal.ssj.rng.*;

/*
 * Tools to prints n points from a two-dimensional point set, for LaTeX figures.
 * Note: There are also several methods to print points in class PointSet.
 * 
 * We also need to be able to add a specific shift mod 1 or a digital shift.
 */
public class PlotPoints {


    /*
     * Print coordinates d1 and d2 of the first n of p, one point per line.
     * Coordinates start at 0.
     */
	public void printPoints (PointSetIterator iter, int d1, int d2, int n) {
        // PointSetIterator iter = p.iterator();
        double[] point = new double[d2+1]; 
        for (int i=0; i < n; i++) {
        	iter.nextPoint (point, d2+1);	
			System.out.println(point[d1] + "  " + point[d2]);
			}
		System.out.println();
		}

	/*
	 * Prints the first two coordinates of all points in p.
	 */
	public void printPoints (PointSet p) {
        PointSetIterator iter = p.iterator();
        printPoints (iter, 0, 1, p.getNumPoints());
    }

	   /*
     * Print points multiplied by n, in binary format
     */
	public void printPointsBinary (PointSetIterator iter, int d1, int d2, int n) {
        double[] point = new double[d2+1]; 
		// System.out.println();
        for (int i=0; i < n; i++) {
        	iter.nextPoint (point, d2+1);    
			System.out.println (PrintfFormat.formatBase (10, 34, 2, point[d1]) 
				+ "   " + PrintfFormat.formatBase (10, 34, 2, point[d2]));
			}
		System.out.println();
		}

  /*
  * Print points multiplied by n, in binary format
  */
	public void printPointsBinaryInt (PointSetIterator iter, int d1, int d2, int n) {
     double[] point = new double[d2+1]; 
		// System.out.println();
     for (int i=0; i < n; i++) {
     	iter.nextPoint (point, d2+1);    
			System.out.println (PrintfFormat.formatBase (8, 2, Math.round(n*point[d1])) 
				+ "   " + PrintfFormat.formatBase (8, 2, Math.round(n*point[d2])));
			}
		System.out.println();
		}

	/*
	 * Prints n pairs of successive pairs from the LCG:  x = ax + c mod m.
	 */
	public void printPointsLcg (long m, long a, long c, long x0, int n) {
		long x = x0;
		long y;
		for (int i=0; i < n; i++) {
			y = (a * x + c) % m;
			System.out.print(x + "  " + y);
			x = y;
			}
		System.out.println();
		}

	
	public static void main(String[] args) {
		
        /*
        * This is a hack to add a fixed digital shift to a DigitalNetBase2.  
        */
		class SingleShiftBase2 extends MRG32k3a {
			// Would need to reset dim, numPoints, etc. ...  !!!!
			int shift0, shift1, shift;   //  31-bit integers

			public void setDigitalShift (int u0, int u1) {
				shift0 = u0;
				shift1 = u1;
				shift = u1;
			}
			public int nextInt (int i, int j) {
				if (shift == shift0) shift = shift1;  
				else shift = shift0;
				return shift;
			}		
		}
		

		PlotPoints plot = new PlotPoints();
        SingleShiftBase2 fixed = new SingleShiftBase2();
        PointSetRandomization rand = new RandomShift(fixed); 
        
        DigitalNetBase2 p = new SobolSequence(6, 31, 2);
        // DigitalNetBase2 p = (new SobolSequence(6, 31, 2)).toNetShiftCj();
        System.out.println (p.toString());
        p.printGeneratorMatrices (2);
        System.out.println ();
        System.out.println (p.formatPoints (p.iteratorNoGray(), 64, 2));

        plot.printPointsBinary(p.iteratorNoGray(), 0, 1, 64);
        plot.printPointsBinary(p.iterator(), 0, 1, 64);
        plot.printPoints(p);
        
        int flip0 = 0b0000000000000000000000000000000;   // Only first bit
        int flip1 = 0b1000000000000000000000000000000;   // Only first bit
        int flip2 = 0b1010000000000000000000000000000;   // Only first bit

        int shift0 = 0b1010010110000111100100110000010;   //  ;
        int shift1 = 0b0101100110001011000100010011010;   //  ;
        // int shift0 = 0b1000000000000000000000000000000;   //  603264427;
        // int shift1 = 0b1000000000000000000000000000000;   //  1074927987;
        
        
        System.out.println ("\n Original binary points, "
        		+ "then with digital shift of 1 bit, then 3 bits.\n");
        
        fixed.setDigitalShift (flip0, flip0);   
        rand.randomize(p);
        plot.printPoints(p.iteratorNoGray(), 0, 1, p.getNumPoints());

        fixed.setDigitalShift (flip1, flip0);   
        rand.randomize(p);
        plot.printPoints(p.iteratorNoGray(), 0, 1, p.getNumPoints());

        fixed.setDigitalShift (flip2, flip0);   
        rand.randomize(p);
        plot.printPoints(p.iteratorNoGray(), 0, 1, p.getNumPoints());

        
        
        fixed.setDigitalShift (shift0, shift1);   
        double dshift0 = (double)shift0 / 2147483648.0;
        double dshift1 = (double)shift1 / 2147483648.0;
		System.out.println (dshift0 + "   " + dshift1);
		System.out.println (shift0 + "   " + shift1);
		System.out.println (PrintfFormat.formatBase (10, 34, 2, dshift0) 
				+ "   " + PrintfFormat.formatBase (10, 34, 2, dshift1));
		System.out.println (PrintfFormat.formatBase (8, 2, shift0) 
				+ "   " + PrintfFormat.formatBase (8, 2, shift1));
		System.out.println();

        
        // PointSetRandomization rand = new LMScrambleShift(new MRG32k3a());
        rand.randomize(p);
        plot.printPointsBinary(p.iteratorNoGray(), 0, 1, p.getNumPoints());
        plot.printPoints(p.iteratorNoGray(), 0, 1, p.getNumPoints());
	}

}
