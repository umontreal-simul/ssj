/*
 * Class:        FaureSequence
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
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
package umontreal.ssj.hups;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class implements digital nets or digital sequences formed by the
 * first @f$n = b^k@f$ points of the Faure sequence in base @f$b@f$. Values
 * of @f$n@f$ up to @f$2^{31}@f$ are allowed. One has @f$r = k@f$. The
 * generator matrices are
 * @f[
 *   \mathbf{C}_j = \mathbf{P}^j \mod b
 * @f]
 * for @f$j=0,…,s-1@f$, where @f$\mathbf{P}@f$ is a @f$k\times k@f$ upper
 * triangular matrix whose entry @f$(l,c)@f$ is the number of combinations of
 * @f$l@f$ objects among @f$c@f$, @f${c\choose l}@f$, for @f$l\le c@f$ and is
 * 0 for @f$l > c@f$. The matrix @f$\mathbf{C}_0@f$ is the identity,
 * @f$\mathbf{C}_1 = \mathbf{P}@f$, and the other @f$\mathbf{C}_j@f$’s can be
 * defined recursively via @f$\mathbf{C}_j = \mathbf{P} \mathbf{C}_{j-1}
 * \mod b@f$. Our implementation uses the recursion
 * @f[
 *   {c \choose l} = {{c-1} \choose l} + {{c-1} \choose{l-1}}
 * @f]
 * to evaluate the binomial coefficients in the matrices @f$\mathbf{C}_j@f$,
 * as suggested by Fox @cite rFOX86a&thinsp; (see also @cite fGLA04a&thinsp;,
 * page 301). The entries @f$x_{j,l,c}@f$ of @f$\mathbf{C}_j@f$ are computed
 * as follows:
 * @f[
 *   \begin{array}{lcll}
 *    x_{j,c,c} 
 *    & 
 *   =
 *    & 
 *    1 
 *    & 
 *   \quad\mbox{ for } c=0,…,k-1,
 *    \\ 
 *    x_{j,0,c} 
 *    & 
 *   =
 *    & 
 *    j x_{j,0,c-1} 
 *    & 
 *   \quad\mbox{ for } c=1,…,k-1, 
 *    \\ 
 *    x_{j,l,c} 
 *    & 
 *   =
 *    & 
 *    x_{j,l-1,c-1} + j x_{j,l,c-1} 
 *    & 
 *   \quad\mbox{ for } 2\le c < l \le k-1, 
 *    \\ 
 *    x_{j,l,c} 
 *    & 
 *   =
 *    & 
 *    0 
 *    & 
 *   \quad\mbox{ for } c>l \mbox{ or } l \ge k. 
 *   \end{array}
 * @f]
 * For any integer @f$m > 0@f$ and @f$\nu\ge0@f$, if we look at the vector
 * @f$(u_{i,j,1},…,u_{i,j,m})@f$ (the first @f$m@f$ digits of coordinate
 * @f$j@f$ of the output) when @f$i@f$ goes from @f$\nu b^m@f$ to
 * @f$(\nu+1)b^m - 1@f$, this vector takes each of its @f$b^m@f$ possible
 * values exactly once. In particular, for @f$\nu= 0@f$, @f$u_{i,j}@f$
 * visits each value in the set @f$\{0, 1/b^m, 2/b^m, …, (b^m-1)/b^m\}@f$
 * exactly once, so all one-dimensional projections of the point set are
 * identical. However, the values are visited in a different order for the
 * different values of @f$j@f$ (otherwise all coordinates would be
 * identical). For @f$j=0@f$, they are visited in the same order as in the
 * van der Corput sequence in base @f$b@f$.
 *
 * An important property of Faure nets is that for any integers @f$m > 0@f$
 * and @f$\nu\ge0@f$, the point set @f$\{\mathbf{u}_i@f$ for @f$i =
 * \nu b^m,…, (\nu+1)b^m -1\}@f$ is a @f$(0,m,s)@f$-net in base @f$b@f$. In
 * particular, for @f$n = b^k@f$, the first @f$n@f$ points form a
 * @f$(0,k,s)@f$-net in base @f$b@f$. The Faure nets are also
 * projection-regular and dimension-stationary (see @cite vLEC02a&thinsp; for
 * definitions of these properties).
 *
 * To obtain digital nets from the *generalized Faure sequence*
 * @cite rTEZ95a&thinsp;, where @f$\mathbf{P}_j@f$ is left-multiplied by some
 * invertible matrix @f$\mathbf{A}_j@f$, it suffices to apply an appropriate
 * matrix scramble (e.g., via  #leftMatrixScramble ). This changes the order
 * in which @f$u_{i,j}@f$ visits its different values, for each coordinate
 * @f$j@f$, but does not change the set of values that are visited. The
 * @f$(0,m,s)@f$-net property stated above remains valid.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class FaureSequence extends DigitalSequence {

    // Maximum dimension for the case where b is not specified.
    // Can be extended by extending the precomputed array prime[].
    private static final int MAXDIM = 500;

    // For storing the generator matrices for given dim and numPoints.
    private int[][][] v;

   /**
    * Constructs a digital net in base @f$b@f$, with @f$n = b^k@f$ points
    * and @f$w@f$ output digits, in `dim` dimensions. The points are the
    * first @f$n@f$ points of the Faure sequence. The generator matrices
    * @f$\mathbf{C}_j@f$ are @f$r\times k@f$. Unless, one plans to apply a
    * randomization on more than @f$k@f$ digits (e.g., a random digital
    * shift for @f$w > k@f$ digits, or a linear scramble yielding @f$r >
    * k@f$ digits), one should take @f$w = r = k@f$ for better
    * computational efficiency. Restrictions: `dim` @f$\le500@f$ and
    * @f$b^k \le2^{31}@f$.
    *  @param b            base
    *  @param k            there will be b^k points
    *  @param r            number of rows in the generator matrices
    *  @param w            number of output digits
    *  @param dim          dimension of the point set
    */
   public FaureSequence (int b, int k, int r, int w, int dim) {
      init (b, k, r, w, dim);
   }

   private void init (int b, int k, int r, int w, int dim) {
      if (dim < 1)
         throw new IllegalArgumentException
            ("Dimension for FaureSequence must be > 1");
      if ((double)k * Math.log ((double)b) > (31.0 * Math.log (2.0)))
         throw new IllegalArgumentException
            ("Trying to construct a FaureSequence with too many points");
      if (r < k || w < r)
         throw new IllegalArgumentException
            ("One must have k <= r <= w for FaureSequence");
      this.b    = b;
      numCols   = k;
      numRows   = r;
      outDigits = w;
      this.dim  = dim;

      int i, j;
      numPoints = b;
      for (i=1; i<k; i++) numPoints *= b;

      // Compute the normalization factors.
      normFactor = 1.0 / Math.pow ((double) b, (double) outDigits);
//      EpsilonHalf = 0.5*normFactor;
      double invb = 1.0 / b;
      factor = new double[outDigits];
      factor[0] = invb;
      for (j = 1; j < outDigits; j++)
         factor[j] = factor[j-1] * invb;

      genMat = new int[dim * numCols][numRows];
      initGenMat();
   }

   /**
    * Same as  {@link #FaureSequence() FaureSequence(b, k, w, w, dim)}
    * with base @f$b@f$ equal to the smallest prime larger or equal to
    * `dim`, and with *at least* `n` points.  The values of @f$k@f$,
    * @f$r@f$, and @f$w@f$ are taken as @f$k = \lceil\log_b n\rceil@f$
    * and @f$r = w = \max(k, \lfloor30 / \log_2 b\rfloor)@f$.
    *  @param n            minimal number of points
    *  @param dim          dimension of the point set
    */
   public FaureSequence (int n, int dim) {
      if ((dim < 1) || (dim > MAXDIM))
         throw new IllegalArgumentException
            ("Dimension for Faure net must be > 1 and < " + MAXDIM);
      b = getSmallestPrime (dim);
      numCols = (int) Math.ceil (Math.log ((double) n)
                                 / Math.log ((double) b));
      outDigits = (int) Math.floor (Math.log ((double)(1 << (MAXBITS - 1)))
                                 / Math.log ((double)b));
      outDigits = Math.max (outDigits, numCols);
      numRows = outDigits;
      init (b, numCols, numRows, outDigits, dim);
   }


   public String toString() {
      StringBuffer sb = new StringBuffer ("Faure sequence:" +
                  PrintfFormat.NEWLINE);
      sb.append (super.toString());
      return sb.toString();
   }


   public void extendSequence (int k) {
      init (b, k, numRows, outDigits, dim);
   }


   // Fills up the generator matrices in genMat for a Faure sequence.
   // See Glasserman (2004), @cite fGLA04a, page 301.
   private void initGenMat() {
      int j, c, l;
      // Initialize C_0 to the identity (for first coordinate).
      for (c = 0; c < numCols; c++) {
         for (l = 0; l < numRows; l++)
            genMat[c][l] = 0;
         genMat[c][c] = 1;
      }
      // Compute C_1, ... ,C_{dim-1}.
      for (j = 1; j < dim; j++) {
         genMat[j*numCols][0] = 1;
         for (c = 1; c < numCols; c++) {
            genMat[j*numCols+c][c] = 1;
            genMat[j*numCols+c][0] = (j * genMat[j*numCols+c-1][0]) % b;
         }
         for (c = 2; c < numCols; c++) {
            for (l = 1; l < c; l++)
               genMat[j*numCols+c][l] = (genMat[j*numCols+c-1][l-1]
                                        + j * genMat[j*numCols+c-1][l]) % b;
         }
         for (c = 0; c < numCols; c++)
            for (l = c+1; l < numRows; l++)
               genMat[j*numCols+c][l] = 0;
      }
/*
      for (j = 0; j < dim; j++) {
     for (l = 0; l < numRows; l++) {
         for (c = 0; c < numCols; c++)
            System.out.print ("  " + genMat[j*numCols+c][l]);
       System.out.println ("");
      }
        System.out.println ("");
  }
*/
   }

/*
   // Fills up the generator matrices in genMat for a Faure net.
   // See Glasserman (2004), @cite fGLA04a, page 301.
   protected void initGenMatNet() {
      int j, c, l, start;
      // Initialize C_0 to the reflected identity (for first coordinate).
      for (c = 0; c < numCols; c++) {
         for (l = 0; l < numRows; l++)
            genMat[c][l] = 0;
         genMat[c][numCols-c-1] = 1;
      }
      // Initialize C_1 to the identity (for second coordinate).
      for (c = 0; c < numCols; c++) {
         for (l = 0; l < numRows; l++)
            genMat[numCols+c][l] = 0;
         genMat[numCols+c][c] = 1;
      }
      // Compute C_2, ... ,C_{dim-1}.
      for (j = 2; j < dim; j++) {
         start = j * numCols;
         genMat[start][0] = 1;
         for (c = 1; c < numCols; c++) {
            genMat[start+c][c] = 1;
            genMat[start+c][0] = ((j-1) * genMat[start+c-1][0]) % b;
         }
         for (c = 2; c < numCols; c++) {
            for (l = 1; l < c; l++)
               genMat[start+c][l] = (genMat[start+c-1][l-1]
                                     + (j-1) * genMat[start+c-1][l]) % b;
         }
         for (c = 0; c < numCols; c++)
            for (l = c+1; l < numRows; l++)
               genMat[start+c][l] = 0;
      }
   }
*/

   // Returns the smallest prime larger or equal to d.
   private int getSmallestPrime (int d) {
      return primes[d-1];
   }

   // Gives the appropriate prime base for each dimension.
   // Perhaps should be internal to getPrime, and non-static, to avoid
   // wasting time and memory when this array is not needed ???
   static final int primes[] =
      {2,2,3,5,5,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,
     23,23,23,29,29,29,29,29,29,31,31,37,37,37,37,37,37,41,41,41,
     41,43,43,47,47,47,47,53,53,53,53,53,53,59,59,59,59,59,59,61,
     61,67,67,67,67,67,67,71,71,71,71,73,73,79,79,79,79,79,79,83,
     83,83,83,89,89,89,89,89,89,97,97,97,97,97,97,97,97,101,101,101,
     101,103,103,107,107,107,107,109,109,113,113,113,113,127,127,127,127,127,127,127,
     127,127,127,127,127,127,127,131,131,131,131,137,137,137,137,137,137,139,139,149,
     149,149,149,149,149,149,149,149,149,151,151,157,157,157,157,157,157,163,163,163,
     163,163,163,167,167,167,167,173,173,173,173,173,173,179,179,179,179,179,179,181,
     181,191,191,191,191,191,191,191,191,191,191,193,193,197,197,197,197,199,199,211,
     211,211,211,211,211,211,211,211,211,211,211,223,223,223,223,223,223,223,223,223,
     223,223,223,227,227,227,227,229,229,233,233,233,233,239,239,239,239,239,239,241,
     241,251,251,251,251,251,251,251,251,251,251,257,257,257,257,257,257,263,263,263,
     263,263,263,269,269,269,269,269,269,271,271,277,277,277,277,277,277,281,281,281,
     281,283,283,293,293,293,293,293,293,293,293,293,293,307,307,307,307,307,307,307,
     307,307,307,307,307,307,307,311,311,311,311,313,313,317,317,317,317,331,331,331,
     331,331,331,331,331,331,331,331,331,331,331,337,337,337,337,337,337,347,347,347,
     347,347,347,347,347,347,347,349,349,353,353,353,353,359,359,359,359,359,359,367,
     367,367,367,367,367,367,367,373,373,373,373,373,373,379,379,379,379,379,379,383,
     383,383,383,389,389,389,389,389,389,397,397,397,397,397,397,397,397,401,401,401,
     401,409,409,409,409,409,409,409,409,419,419,419,419,419,419,419,419,419,419,421,
     421,431,431,431,431,431,431,431,431,431,431,433,433,439,439,439,439,439,439,443,
     443,443,443,449,449,449,449,449,449,457,457,457,457,457,457,457,457,461,461,461,
     461,463,463,467,467,467,467,479,479,479,479,479,479,479,479,479,479,479,479,487,
     487,487,487,487,487,487,487,491,491,491,491,499,499,499,499,499,499,499,499,503};

}