/*
 * Class:        HilbertCurveMap
 * Description:  Map the Hilbert curve in a d-dimensional space [0,1)^d.
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2014  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since

 * SSJ is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, either version 3 of the License, or
 * any later version.

 * SSJ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * A copy of the GNU General Public License is available at
   <a href="http://www.gnu.org/licenses">GPL licence site</a>.
 */

  /* IMPORTANT NOTE:
	* Much of this code has been taken (with adaptations) from  
  *     the hilbert.c  code  
  * Author:	Spencer W. Thomas
  * 		EECS Dept.
  * 		University of Michigan
  * Date:	Thu Feb  7 1991
  * Copyright (c) 1991, University of Michigan
  */
package umontreal.ssj.util.sort;
  import java.util.Comparator;
  import java.util.Arrays;

/**
 * This class implements the mapping of a Hilbert curve in the
 * @f$d@f$-dimensional unit hypercube @f$[0,1)^d@f$. This mapping can be used
 * for sorting algorithms.
 *
 * This map (conceptually) divides the unit hypercube @f$[0,1)^d@f$ in
 * @f$2^{dm}@f$ subcubes of equal sizes, by dividing each axis in @f$2^m@f$
 * equal parts, and uses the first @f$m@f$ bits of each of the @f$d@f$
 * coordinates to place each point in one of the subcubes. It then enumerates
 * the subcubes in the same order as a Hilbert curve in @f$[0,1)^d@f$ would
 * visit them, and orders the points accordingly. Each cube has an (integer)
 * Hilbert index @f$r@f$ from 0 to @f$2^{dm}-1@f$ and the cubes (and points)
 * are ordered according to this index. Two points that fall in the same
 * subcube can be placed in an unspecified (arbitrary) order.
 *
 * For the implementation of sorts based on the Hilbert curve (or Hilbert
 * index), we identify and sort the subcubes by their Hilbert index, but it
 * is also convenient to identify them (alternatively) with @f$m@f$-bit
 * integer coordinates: The subcube with coordinates @f$(i_1,…,i_d)@f$ is
 * defined as @f$\prod_{j=0}^{d-1} [i_j 2^{-m},  (i_j+1) 2^{-m})@f$. Note
 * that each interval is open on the right. That is, if we multiply the
 * coordinates of a point in the subcube by @f$2^m@f$ and truncate them to
 * integers, we obtain the integer coordinates of the subcube. For example,
 * if @f$d=2@f$ and @f$m=4@f$, we have @f$2^8 = 256@f$ subcubes, whose
 * integer coordinates go from 0 to 15, and the point @f$(0.1, 0.51)@f$
 * belongs to the subcube with integer coordinates @f$(1, 8)@f$.
 *
 * For given @f$d@f$ and @f$m@f$, this class offers methods to compute the
 * integer coordinates of the corresponding subcube from the real-valued
 * coordinates of a point in @f$[0,1)^d@f$, as well as the Hilbert index of a
 * subcube from its integer coordinates, and vice versa. The code that
 * computes the latter correspondences is taken (with slight adaptations)
 * from the `hilbert.c` program of Spencer W. Thomas, University of Michigan,
 * 1991.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class HilbertCurveMap {

    int dimension;  // Dimension d of the points used for the sort. 
    int m;       // Number of bit retained for each coordinate. Must be <= 31.
    // int twom;    // Number of intervals per coordinate = 2^{m}.   
    // long nCubes;   // Total number of subcubes = 2^{dm}.   Needed?   **********
    int maxnbits;    // Needed ???  It seems this should be = dimension    ******
    int maxlength;   // 2^d  =  (1 << d)

    /*
 	* The following are precomputed tables to simplify calculations.  
	* Notation: p#i means bit i in byte p (high order bit first).
 	* p_to_s:	Output s is a byte from input p such that
 	* 		s#i = p#i xor p#(i-1)
 	* s_to_p:	The inverse of the above.
 	* p_to_J:	Output J is "principle position" of input p.  The
 	* 		principle position is the last bit s.t.
 	* 		p#J != p#(n-1) (or n-1 if all bits are equal).
 	* bit:		bit[i] == (1 << (n - i))
 	* circshift:	circshift[b][i] is a right circular shift of b by i
 	* 		bits in n bits.
 	* parity:	Parity[i] is 1 or 0 depending on the parity of i (1 is odd).
 	* bitof:	bitof[b][i] is b#i.
 	* nbits:	The value of n for which the above tables have been
 	* 		calculated.
 	*/


    // int dimForTables = 0; // Dimensions in which tables have been computed.
    int[] p_to_s;
    int[] s_to_p;
    int[] p_to_J;
    int[] bit;
    int[] parity;
    int[][] circshift;
    int[][] bitof;


    // Precomputes several tables.  
    // This code is from Spencer W. Thomas.
    void initTables () {
        p_to_s = new int[maxlength];
        s_to_p = new int[maxlength];
        p_to_J = new int[maxlength];
        parity = new int[maxlength];
        bit    = new int[maxnbits];
        circshift = new int[maxlength][maxnbits];
        bitof     = new int[maxlength][maxnbits];
        //  It seems that maxnbits  should be the dimension d !!!

        int i, b;
				int n = dimension;   // We use n to avoid changing this code.
        long two_n = 1 << n;

        // if (dimForTables == n )
        //    return;       // Tables have been computed already!
        // dimForTables = n;
        /* bit array is easy. */
        for ( b = 0; b < n; b++ )
            bit[b] = (int)(1 << (n - b - 1));

        /* Next, do bitof. */
        for ( i = 0; i < two_n; i++ )
            for ( b = 0; b < n; b++ )
                bitof[i][b] = (int) ((i & bit[b]) !=0 ? 1 : 0);

        /* circshift is independent of the others. */
        for ( i = 0; i < two_n; i++ )
            for ( b = 0; b < n; b++ )
                circshift[i][b] = (int) ((i >> (b)) |
                        ((i << (n - b)) & (two_n - 1)));

        /* So is parity. */
        parity[0] = 0;
        for ( i = 1, b = 1; i < two_n; i++ )
        {
            /* Parity of i is opposite of the number you get when you
             * knock the high order bit off of i.
             */
            if ( i == b * 2 )
                b *= 2;
            parity[i] = (int) (1 - parity[i - b]);
        }

        /* Now do p_to_s, s_to_p, and p_to_J. */
        for ( i = 0; i < two_n; i++ )
        {
            int s;
            s = i & bit[0];
            for ( b = 1; b < n; b++ )
                if ( (bitof[i][b] ^ bitof[i][b-1]) != 0)
                    s |= bit[b];
            p_to_s[i] = (int) s;
            s_to_p[s] = (int) i;

            p_to_J[i] = (int) (n - 1);
            for ( b = 0; b < n; b++ )
                if ( bitof[i][b] != bitof[i][n-1] )
                    p_to_J[i] = (int) b;
        }
    }

   /**
    * Constructs a  @ref HilbertCurveMap object that will use the first
    * @f$m@f$ bits of each of the first `d` coordinates to sort the
    * points. In this implementation, one must have @f$md \le63@f$.
    *  @param d            maximum dimension
    *  @param m            number of bits used for each coordinate
    */
   public HilbertCurveMap (int d, int m) {      
      dimension = d;
      this.m = m;
      // twom = 1 << m;
      // nCubes = 1 << d * m;

      maxnbits = m;    //    ???            ****************
      maxlength = 1 << dimension;    // 2^d
      this.initTables();
    }

   /**
    * This constructor is similar to  #HilbertCurveMap(int,int), but the
    * parameter @f$m = \lfloor63/d\rfloor@f$. In this implementation,
    * one must have @f$md \le63@f$.
    *  @param d            maximum dimension
    */
   public HilbertCurveMap (int d) {      
      this (d, 63 / d);
    }

    /**
     * Returns the dimension of the unit hypercube.
     */
    public int dimension() {
       return dimension;
    }

    /**
     * Returns the number of bits @f$m@f$ that is used to divide the axis
     * of each coordinate.
     */
    public int getM() {
       return m;
    }


    // The following code is adapted from hilbert.c, by Spencer W. Thomas.

    /*****************************************************************
     * TAG( hilbert_i2c )
     * 
     * Convert an index into a Hilbert curve to a set of coordinates.
     * Inputs:
     * 	n:	Number of coordinate axes.   (= d)
     * 	m:	Number of bits per axis.
     * 	r:	The index, contains n*m bits (so n*m must be <= 63).
     * Outputs:
     * 	a:	The list of n coordinates, each with m bits.
     * Assumptions:
     * 	n*m < (sizeof r) * (bits_per_int).
     * Algorithm:
     * 	From A. R. Butz, "Alternative Algorithm for Hilbert's
     * 		Space-Filling Curve", IEEE Trans. Comp., April, 1971,
     * 		pp 424-426.
     */

/**
 * Takes as input the Hilbert index @f$r@f$ of a subcube and returns in `a[]`
 * its integer coordinates. WARNING: This method currently works only for
 * @f$m \le9@f$. It is not used for the sort.
 */
public void indexToCoordinates (long r, int[] a) {
        int[] rho = new int[9];    //  What is this 9 doing here ???      *******
        int rh, J, sigma, tau;
        int sigmaT, tauT; 
        int tauT1 = 0;
        int omega, omega1 = 0;
        int[] alpha = new int[9];  //  What is this 9 doing here ???      *******
        int i, b;
        int Jsum;
				int n = dimension;

        /* Distribute bits from r into rho. */
        for ( i = m - 1; i >= 0; i-- )
        {
            rho[i] = (int) (r & ((1 << n) - 1));
            r >>= n;
        }

        /* Loop over ints. */
        Jsum = 0;
        for ( i = 0; i < m; i++ )
        {
            rh = rho[i];
            /* J[i] is principle position of rho[i]. */
            J = p_to_J[rh];

            /* sigma[i] is derived from rho[i] by exclusive-oring adjacent bits. */
            sigma = p_to_s[rh];

            /* tau[i] complements low bit of sigma[i], and bit at J[i] if
             * necessary to make even parity.
             */
            tau = (int) (sigma ^ 1);
            if ( parity[tau] != 0 )
                tau ^= bit[J];

            /* sigmaT[i] is circular shift of sigma[i] by sum of J[0..i-1] */
            /* tauT[i] is same circular shift of tau[i]. */
            if ( Jsum > 0 )
            {
                sigmaT = circshift[sigma][Jsum];
                tauT = circshift[tau][Jsum];
            }
            else
            {
                sigmaT = sigma;
                tauT = tau;
            }

            Jsum += J;
            if ( Jsum >= n )
                Jsum -= n;

            /* omega[i] is xor of omega[i-1] and tauT[i-1]. */
            if ( i == 0 )
                omega = 0;
            else
                omega = (int) (omega1 ^ tauT1);
            omega1 = omega;
            tauT1 = tauT;

            /* alpha[i] is xor of omega[i] and sigmaT[i] */
            alpha[i] = (int) (omega ^ sigmaT);
        }

        /* Build coordinates by taking bits from alphas. */
        for ( b = 0; b < n; b++ )
        {
            int ab, bt;
            ab = 0;
            bt = bit[b];
            /* Unroll the loop that stuffs bits into ab.
             * The result is shifted left by 9-m bits.
             */
            
            switch( m )
            {
                case 9:	if ( (alpha[8] & bt) != 0) ab |= 0x01;
                case 8:	if ( (alpha[7] & bt) != 0) ab |= 0x02;
                case 7:	if ( (alpha[6] & bt) != 0) ab |= 0x04;
                case 6:	if ( (alpha[5] & bt) != 0) ab |= 0x08;
                case 5:	if ( (alpha[4] & bt) != 0) ab |= 0x10;
                case 4:	if ( (alpha[3] & bt) != 0) ab |= 0x20;
                case 3:	if ( (alpha[2] & bt) != 0) ab |= 0x40;
                case 2:	if ( (alpha[1] & bt) != 0) ab |= 0x80;
                case 1:	if ( (alpha[0] & bt) != 0) ab |= 0x100;
            }
            a[b] = ab >> (9 - m);
        }
    }


        /*****************************************************************
     * TAG( hilbert_c2i )
     * 
     * Convert coordinates of a point on a Hilbert curve to its index.
     * Inputs:
     * 	n:	Number of coordinates.    (= d)
     * 	m:	Number of bits/coordinate.
     * 	a:	Array of n m-bit coordinates.
     * Outputs:
     * 	r:	Output index value.  n*m bits.
     * Assumptions:
     * 	n*m <= 63.
     * Algorithm:
     * 	Invert the above.
     */

/**
 * Takes as input in `a[]` the integer coordinates of a subcube and returns
 * the corresponding Hilbert index @f$r@f$.
 */
public long coordinatesToIndex (int a[]) {
        int[] rho = new int[maxnbits];
        int[] alpha = new int[maxnbits];
        int J, sigma, tau, sigmaT, tauT, tauT1 = 0, omega, omega1 = 0;
        int i, b;
        int Jsum;
        long rl;
        int n = dimension;
								
        //  initTables();
 
        /* Unpack the coordinates into alpha[i]. */
        /* First, zero out the alphas. */
        for(i = m; i>0; --i)
        {
            alpha[i-1] = 0;
        }
        for ( b = 0; b < n; b++ )
        {
            long bt = bit[b];
            long t = a[b];
            for(i = 1; i <= m; ++i)
            {
                if ((t >> (m-i)) != 0) 
                {
                    alpha[i-1] |= bt;
                    t = t - (1 <<(m-i)); 
                }
            }
        }

        Jsum = 0;
        for ( i = 0; i < m; i++ )
        {
            /* Compute omega[i] = omega[i-1] xor tauT[i-1]. */
            if ( i == 0 )
                omega = 0;
            else
                omega = (int) (omega1 ^ tauT1);

            sigmaT = (int) (alpha[i] ^ omega);
            /* sigma[i] is the left circular shift of sigmaT[i]. */
            if ( Jsum != 0 )
                sigma = circshift[sigmaT][n - Jsum];
            else
                sigma = sigmaT;

            rho[i] = s_to_p[sigma];

            /* Now we can get the principle position. */
            J = p_to_J[rho[i]];

            /* And compute tau[i] and tauT[i]. */
            /* tau[i] complements low bit of sigma[i], and bit at J[i] if
             * necessary to make even parity.
             */
            tau = (int) (sigma ^ 1);
            if ( parity[tau] !=0 )
                tau ^= bit[J];

            /* tauT[i] is right circular shift of tau[i]. */
            if ( Jsum != 0 )
                tauT = circshift[tau][Jsum];
            else
                tauT = tau;
            Jsum += J;
            if ( Jsum >= n )
                Jsum -= n;

            /* Carry forth the "i-1" values. */
            tauT1 = tauT;
            omega1 = omega;
        }

        /* Pack rho values into r. */
        rl = 0;
        for ( i = 0; i < m; i++ )
            rl = (rl << n) | rho[i];
        return rl;
  }

    /**
     * Takes in `p[]` a point with real-valued coordinates and places in
     * `a[]` the integer coordinates of the corresponding subcube. The
     * two arrays are assumed to have length @f$d@f$.
     */
    public void pointToCoordinates (double[] p, int[] a) {
       for (int j = 0; j < dimension; ++j)
          a[j] = (int) (p[j] * (1 << m));  // Multiply by 2^m.      
    }

}
