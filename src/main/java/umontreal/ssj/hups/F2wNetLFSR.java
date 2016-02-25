/*
 * Class:        F2wNetLFSR
 * Description:  digital net in base 2 starting from a linear feedback
                 shift register generator
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
 * This class implements a digital net in base 2 starting from a linear
 * feedback shift register generator. It is exactly the same point set as the
 * one defined in the class  @ref umontreal.ssj.hups.F2wCycleBasedLFSR. See
 * the description of this class for more details on the way the point set is
 * constructed.
 *
 * Constructing a point set using this class, instead of using
 * @ref umontreal.ssj.hups.F2wCycleBasedLFSR, makes SSJ construct a digital
 * net in base 2. This is useful if one wants to randomize using one of the
 * randomizations included in the class  @ref umontreal.ssj.hups.DigitalNet.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class F2wNetLFSR extends DigitalNetBase2 
{
   private F2wStructure param;

    /**
     * Constructs and stores the set of cycles for an LCG with
     *    modulus <SPAN CLASS="MATH"><I>n</I></SPAN> and multiplier <SPAN CLASS="MATH"><I>a</I></SPAN>.
     *   If pgcd<SPAN CLASS="MATH">(<I>a</I>, <I>n</I>) = 1</SPAN>, this constructs a full-period LCG which has two
     *   cycles, one containing 0 and one, the LCG period.
     *
     * @param n required number of points and modulo of the LCG
     *
     *    @param a generator <SPAN CLASS="MATH"><I>a</I></SPAN> of the LCG
     *
     *
     */

   /**
    * Constructs a point set with @f$2^{rw}@f$ points. See the description
    * of the class  @ref umontreal.ssj.hups.F2wStructure for the meaning
    * of the parameters.
    */
   public F2wNetLFSR (int w, int r, int modQ, int step, int nbcoeff,
                      int coeff[], int nocoeff[], int dim) 
   {
      param = new F2wStructure (w, r, modQ, step, nbcoeff, coeff, nocoeff);
      param.initParamLFSR ();
      initNet (r, w, dim);
   }

   /**
    * Constructs a point set after reading its parameters from file
    * `filename`; the parameters are located at line numbered `no` of
    * `filename`. The available files are listed in the description of
    * class  @ref umontreal.ssj.hups.F2wStructure.
    */
   public F2wNetLFSR (String filename, int no, int dim) 
   {
      param = new F2wStructure (filename, no);
      param.initParamLFSR ();
      initNet (param.r, param.w, dim);
   }


   public String toString ()
   {
       String s = "F2wNetLFSR:" + PrintfFormat.NEWLINE;
       return s + param.toString ();
   }


   private void initNet (int r, int w, int dim)
   {
      numCols = r * w;
      numRows = 31;
      outDigits = 31;
      numPoints = (1 << numCols);
      this.dim = dim;
      normFactor = 1.0 / (1L << 31);
      genMat = new int[dim * numCols];

      for (int j = 0; j < numCols; j++) {
         param.state = 1 << (r * w - 1 - j);
         param.initF2wLFSR ();
         genMat[j] = param.output;
         for (int i = 1; i < dim; i++) {
            param.F2wLFSR ();
            genMat[i * numCols + j] = param.output;
         }
      }
   }
}