/*
 * Class:        NiedSequenceBase2
 * Description:  digital Niederreiter sequences in base 2.
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
 * This class implements digital sequences constructed from the Niederreiter
 * sequence in base 2. For details on these point sets, see
 * @cite rBRA92a&thinsp;.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class NiedSequenceBase2 extends DigitalSequenceBase2 { 

   private static final int MAXDIM = 318;  // Maximum dimension.
   private static final int NUMCOLS = 30;  // Maximum number of columns.

   /**
    * Constructs a new digital sequence in base 2 from the first
    * @f$n=2^k@f$ points of the Niederreiter sequence, with @f$w@f$ output
    * digits, in `dim` dimensions. The generator matrices
    * @f$\mathbf{C}_j@f$ are @f$w\times k@f$. Restrictions:
    * @f$0\le k\le30@f$, @f$k\le w@f$, and `dim` @f$\le318@f$.
    *  @param k            there will be 2^k points
    *  @param w            number of output digits
    *  @param dim          dimension of the point set
    */
   public NiedSequenceBase2 (int k, int w, int dim) {
      init (k, w, w, dim);
   }


   public String toString() {
      StringBuffer sb = new StringBuffer ("Niederreiter sequence in base 2" +
                                           PrintfFormat.NEWLINE);
      sb.append (super.toString());
      return sb.toString();
   }

   private void init (int k, int r, int w, int dim) {
      if ((dim < 1) || (dim > MAXDIM))
         throw new IllegalArgumentException 
            ("Dimension for NiedSequenceBase2 must be > 1 and <= " + MAXDIM);
      if (r < k || w < r || w > MAXBITS || k >= MAXBITS) 
         throw new IllegalArgumentException
            ("One must have k < 31 and k <= r <= w <= 31 for NiedSequenceBase2");
      numCols   = k;
      numRows   = r;   // Unused!
      outDigits = w;
      numPoints = (1 << k);
      this.dim  = dim;
      normFactor = 1.0 / ((double) (1L << (outDigits)));
      genMat = new int[dim * numCols];
      initGenMat();
   }


   public void extendSequence (int k) {
      init (k, numRows, outDigits, dim);
   }


   // Initializes the generator matrices for a sequence. 
   /* I multiply by 2 because the relevant columns are in the 30 least 
      significant bits of NiedMat, but if I understand correctly,
      SSJ assumes that they are in bits [31, ..., 1]. 
      Then I shift right if w < 31. */
   private void initGenMat ()  {
      for (int j = 0; j < dim; j++)
         for (int c = 0; c < numCols; c++) {
            genMat[j*numCols + c] = NiedMat[j*NUMCOLS + c] << 1;
            genMat[j*numCols + c] >>= MAXBITS - outDigits;
         }
   }

/*
   // Initializes the generator matrices for a net. 
   protected void initGenMatNet()  {
      int j, c;

      // the first dimension, j = 0.
      for (c = 0; c < numCols; c++)
         genMat[c] = (1 << (outDigits-numCols+c));

      for (j = 1; j < dim; j++)
         for (c = 0; c < numCols; c++)
            genMat[j*numCols + c] = 2 * NiedMat[(j - 1)*NUMCOLS + c];
   }
*/

   // ****************************************************************** 
   // Generator matrices of Niederreiter sequence. 
   // This array stores explicitly NUMCOLS columns in 318 dimensions.

   private static int[] NiedMat;
   private static final int MAXN = 9540;

   static {
      NiedMat = new int[MAXN];

      try {
         InputStream is = 
            NiedSequenceBase2.class.getClassLoader().getResourceAsStream (
            "umontreal/ssj/hups/dataSer/Nieder/NiedSequenceBase2.ser");
         if (is == null)
            throw new FileNotFoundException (
               "Cannot find NiedSequenceBase2.ser");
         ObjectInputStream ois = new ObjectInputStream(is);
         NiedMat = (int[]) ois.readObject();
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