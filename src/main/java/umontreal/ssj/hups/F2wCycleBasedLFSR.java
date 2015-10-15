/*
 * Class:        F2wCycleBasedLFSR
 * Description:  point set based upon a linear feedback shift register
                 sequence
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
import cern.colt.list.*;
import umontreal.ssj.util.PrintfFormat;

/**
 * This class creates a point set based upon a linear feedback shift register
 * sequence. The recurrence used to produce the point set is
 * @f[
 *   m_n = \sum_{i=1}^r b_i m_{n-i}
 * @f]
 * where @f$m_n\in\mathbb F_{2^w}@f$, @f$n\geq0@f$ and @f$b_i\in\mathbb
 * F_{2^w}@f$. There is a polynomial in @f$\mathbb F_{2^w}[z]@f$ associated
 * with this recurrence called the *characteristic polynomial*. It is
 * @f[
 *   P(z) = z^r + \sum_{i=1}^r b_i z^{r-i}.
 * @f]
 * In the implementation, this polynomial is stored in an object
 * `F2wStructure`.
 *
 * Let @f${\mathbf{x}} = (x^{(0)}, …, x^{(p-1)}) \in\mathbb F_2^p@f$ be a
 * @f$p@f$-bit vector. Let us define the function @f$\phi(\mathbf{x}) =
 * \sum_{i=1}^p 2^{-i} x^{(i-1)}@f$. The point set in @f$t@f$ dimensions
 * produced by this class is
 * @f[
 *   \left\{ (\phi(\mathbf{y}_0),\phi(\mathbf{y}_s),…,\phi(\mathbf{y}_{s(t-1)}): (\mathbf{v}_0,…,\mathbf{v}_{r-1})\in\mathbb F_2^{rw}\right\}
 * @f]
 * where @f$\mathbf{y}_n = \mbox{trunc}_h(\mathbf{v}_n,
 * \mathbf{v}_{n+1},…)@f$, @f$\mathbf{v}_n@f$ is the representation of
 * @f$m_n@f$ under the polynomial basis of @f$\mathbb F_{2^w}@f$ over
 * @f$\mathbb F_2@f$, and @f$h=w\lfloor31/w\rfloor@f$. The parameter
 * @f$s@f$ is called the stepping parameter of the recurrence.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class F2wCycleBasedLFSR extends CycleBasedPointSetBase2 {

   private F2wStructure param;

   /**
    * Constructs a point set with @f$2^{rw}@f$ points. See the description
    * of the class  @ref umontreal.ssj.hups.F2wStructure for the meaning
    * of the parameters.
    */
   public F2wCycleBasedLFSR (int w, int r, int modQ, int step, int nbcoeff,
                             int coeff[], int nocoeff[])
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
   {
      param = new F2wStructure (w, r, modQ, step, nbcoeff, coeff, nocoeff);
      init ();
   }

   /**
    * Constructs a point set after reading its parameters from file
    * `filename`; the parameters are located at line numbered `no` of
    * `filename`. The available files are listed in the description of
    * class  @ref umontreal.ssj.hups.F2wStructure.
    */
   public F2wCycleBasedLFSR (String filename, int no)
   {
      param = new F2wStructure (filename, no);
      init ();
   }


   private void init ()
   {
      param.initParamLFSR ();
      normFactor = param.normFactor;
      EpsilonHalf = param.EpsilonHalf;
      numBits = param.numBits;
      fillCyclesLFSR ();
   }

   public String toString ()
   {
       String s = "F2wCycleBasedLFSR:" + PrintfFormat.NEWLINE;
       return s + param.toString ();
   }

   private void fillCyclesLFSR ()
   {
      int n = 1 << param.getLog2N ();
      IntArrayList c;             // Array used to store the current cycle.
      int currentState;           // The state currently visited.
      int i;
      boolean stateVisited[] = new boolean[n];

      // Indicates which states have been visited so far.
      for (i = 0; i < n; i++)
         stateVisited[i] = false;
      int startState = 0;    // First state of the cycle currently considered.
      numPoints = 0;
      while (startState < n) {
         stateVisited[startState] = true;
         c = new IntArrayList ();
         param.state = startState;
         param.initF2wLFSR ();
         c.add (param.output);
         param.F2wLFSR ();
         // Warning: watch for overflow !!!
         while (param.state != startState) {
            stateVisited[param.state] = true;
            c.add (param.output);
            param.F2wLFSR ();
         }
         addCycle (c);
         for (i = startState + 1; i < n; i++)
            if (stateVisited[i] == false)
               break;
         startState = i;
      }
   }
}