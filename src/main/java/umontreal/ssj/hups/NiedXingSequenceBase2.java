/*
 * Class:        NiedXingSequenceBase2
 * Description:  Niederreiter-Xing sequences in base 2
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

import java.io.Serializable;
import java.io.ObjectInputStream;
import java.io.InputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class implements digital sequences based on the Niederreiter-Xing
 * sequence in base 2. For details on these point sets, see @cite rNIE98a,
 * @cite rPIR01c&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class NiedXingSequenceBase2 extends DigitalSequenceBase2 { 

   private static final int MAXDIM  = 32;  // Maximum dimension.
   private static final int NUMCOLS = 30;  // Maximum number of columns.
   private static final boolean isTrans = true;

   /**
    * Constructs a new Niederreiter-Xing digital sequence in base 2 with
    * @f$n = 2^k@f$ points and @f$w@f$ output digits, in `dim` dimensions.
    * The generator matrices @f$\mathbf{C}_j@f$ are @f$w\times k@f$ and
    * the numbers making the bit matrices are taken from [Pirsic’s
    * site](http://www.ricam.oeaw.ac.at/people/page/pirsic/niedxing/index.html)
    * . The bit matrices from Pirsic’s site are transposed to be
    * consistent with SSJ, and at most 30 bits of the matrices are used.
    * Restrictions: @f$0\le k\le30@f$, @f$k\le w@f$, and `dim`
    * @f$\le32@f$.
    *  @param k            there will be 2^k points
    *  @param w            number of output digits
    *  @param dim          dimension of the point set
    */
   public NiedXingSequenceBase2 (int k, int w, int dim) {
      init (k, w, w, dim);
   }


   public String toString() {
      StringBuffer sb = new StringBuffer ("Niederreiter-Xing sequence:" +
                                           PrintfFormat.NEWLINE);
      sb.append (super.toString());
      return sb.toString();
   }

   private void init (int k, int r, int w, int dim) {
      if ((dim < 1) || (dim > MAXDIM))
         throw new IllegalArgumentException 
            ("Dimension for NiedXingSequenceBase2 must be > 1 and <= " + MAXDIM);
      if (r < k || w < r || w > MAXBITS || k >= MAXBITS) 
         throw new IllegalArgumentException
        ("One must have k < 31 and k <= r <= w <= 31 for NiedXingSequenceBase2");
      numCols   = k;
      numRows   = r;
      outDigits = w;
      numPoints = (1 << k);
      this.dim  = dim;
      // 1L otherwise gives wrong results for outDigits >= 31
      normFactor = 1.0 / ((double) (1L << (outDigits)));
      genMat = new int[dim * numCols];
      initGenMat();
   } 


   public void extendSequence (int k) {
      init (k, numRows, outDigits, dim);
   }


   // Initializes the generator matrices for a sequence.
   private void initGenMat ()  {
      // Compute where we should start reading matrices.
      /* Pirsic gives matrices starting at dimension 4; there are j matrices
         of dimension j. Thus if we want to use matrices in dimension
         dim, we must jump over all matrices of dimension < dim. If they
         started at dimension 1, there would be dim*(dim-1)/2 elements to
         disregard. But the first 3 dimensions are not there, thus subtract
         3*(3 - 1)/2 = 6 to get the starting element of dimension dim.
      I also multiply by 2 because the relevant columns are in the 30 least 
         significant bits of NiedXingMat, but if I understand correctly,
         SSJ assumes that they are in bits [31, ..., 1].
      At last, I shift right if w < 31 bits. */

      int start;
      if (dim <= 4) 
         start = 0;
      else
         start = ((dim*(dim-1)/2)-6)*NUMCOLS;

      long x;
      if (isTrans) {
         for (int j = 0; j < dim; j++)
            for (int c = 0; c < numCols; c++) {
               x = NiedXingMatTrans[start + j*NUMCOLS + c];
               x <<= 1;
               genMat[j*numCols + c] = (int) (x >> (MAXBITS - outDigits));
             }
      } else {
         for (int j = 0; j < dim; j++)
            for (int c = 0; c < numCols; c++) {
               x = NiedXingMat[start + j*NUMCOLS + c];
               x <<= 1;
               genMat[j*numCols + c] = (int) (x >> (MAXBITS - outDigits));
            }
      }
   }


   // ****************************************************************** 
   // Generator matrices of Niederreite-Xing sequences.
   // This array stores explicitly NUMCOLS columns in MAXDIM dimensions.
   // The implemented generator matrices are provided by Gottlieb Pirsic
   // and were downloaded on 15 July 2004. (RS)
   // These are the numbers given by Pirsic: I kept the first 30 of each
   // row vector in each dimension and shifted them to get the 30 most 
   // significant bits. These numbers considered as 30 X 30 bit matrices 
   // are given in matrix NiedXingMat below. The same numbers but transposed
   // as 30 X 30 bit matrices are given in matrix NiedXingMatTrans below.
   // According to Yves Edel, the correct matrices for Niederreiter-Xing
   // are NiedXingMatTrans.
   //
   // The matrices were given up to MAXDIM = 32, but the javac compiler
   // cannot compile code  that is too big. So I serialized them in
   // file NiedXingSequenceBase2.ser. (RS)

   private static int[] NiedXingMat;
   private static int[] NiedXingMatTrans;
   private static final int MAXN = 15660;

   static {
//      NiedXingMat = new int[MAXN];
      NiedXingMatTrans = new int[MAXN];

      try {
/*
         InputStream is = 
            NiedXingSequenceBase2.class.getClassLoader().getResourceAsStream (
         "umontreal/ssj/hups/dataSer/Nieder/NiedXingSequenceBase2.ser");
         if (is == null)
            throw new FileNotFoundException (
               "Cannot find NiedXingSequenceBase2.ser");
         ObjectInputStream ois = new ObjectInputStream(is);
         NiedXingMat = (int[]) ois.readObject();
         ois.close();
*/
         InputStream is =
            NiedXingSequenceBase2.class.getClassLoader().getResourceAsStream(
      "umontreal/ssj/hups/dataSer/Nieder/NiedXingSequenceBase2Trans.ser");
         if (is == null)
            throw new FileNotFoundException (
               "Cannot find NiedXingSequenceBase2Trans.ser");
         ObjectInputStream ois = new ObjectInputStream(is);
         NiedXingMatTrans = (int[]) ois.readObject();
         ois.close();

      } catch(FileNotFoundException e) {
         e.printStackTrace();
         System.exit(1);

      } catch(IOException e) {
         e.printStackTrace();
         System.exit(1);

      } catch(ClassNotFoundException e) {
         e.printStackTrace();
         System.exit(1);
      }
   }
}