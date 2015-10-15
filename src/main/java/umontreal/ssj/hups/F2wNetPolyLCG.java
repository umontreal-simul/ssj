/*
 * Class:        F2wNetPolyLCG
 * Description:  digital nets in base 2 starting from a polynomial LCG 
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
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
package umontreal.ssj.hups; 

import umontreal.ssj.util.PrintfFormat;

/**
 * This class implements a digital net in base 2 starting from a polynomial
 * LCG in @f$\mathbb F_{2^w}[z]/P(z)@f$. It is exactly the same point set as
 * the one defined in the class
 * @ref umontreal.ssj.hups.F2wCycleBasedPolyLCG. See the description of this
 * class for more details on the way the point set is constructed.
 *
 * Constructing a point set using this class, instead of using
 * @ref umontreal.ssj.hups.F2wCycleBasedPolyLCG, makes SSJ construct a
 * digital net in base 2. This is useful if one wants to randomize using one
 * of the randomizations included in the class
 * @ref umontreal.ssj.hups.DigitalNet.
 *
 * **Note: This class in not operational yet!**
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class F2wNetPolyLCG extends DigitalNetBase2 
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
   public F2wNetPolyLCG (int type, int w, int r, int modQ, int step,
                         int nbcoeff, int coeff[], int nocoeff[], int dim) 
   {
      param = new F2wStructure (w, r, modQ, step, nbcoeff, coeff, nocoeff);
      initNet (r, w, dim);
   }

   /**
    * Constructs a point set after reading its parameters from file
    * `filename`; the parameters are located at line numbered `no` of
    * `filename`. The available files are listed in the description of
    * class  @ref umontreal.ssj.hups.F2wStructure.
    */
   public F2wNetPolyLCG (String filename, int no, int dim) 
   {
      param = new F2wStructure (filename, no);
      initNet (param.r, param.w, dim);
   }
 

   public String toString ()
   {
       String s = "F2wNetPolyLCG:" + PrintfFormat.NEWLINE;
       return s + param.toString ();
   }


   private void initNet (int r, int w, int dim)
   {
      normFactor = param.normFactor;
   }
}