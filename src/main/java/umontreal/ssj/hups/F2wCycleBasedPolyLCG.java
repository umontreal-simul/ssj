/*
 * Class:        F2wCycleBasedPolyLCG
 * Description:  point set based upon a linear congruential sequence in a
                 finite field
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
import umontreal.ssj.rng.*;
import cern.colt.list.*;

/**
 * This class creates a point set based upon a linear congruential sequence
 * in the finite field @f$\mathbb F_{2^w}[z]/P(z)@f$. The recurrence is
 * @f[
 *   q_n(z) = z^s q_{n-1}(z) \mod P(z)
 * @f]
 * where @f$P(z)\in\mathbb F_{2^w}[z]@f$ has degree @f$r@f$ and @f$q_n(z) =
 * q_{n,1} z^{r-1} + \cdots+ q_{n,r} \in\mathbb F_{2^w}[z]/P(z)@f$. The
 * parameter @f$s@f$ is called the stepping parameter of the recurrence. The
 * polynomial @f$P(z)@f$ is not necessarily the characteristic polynomial of
 * this recurrence, but it can still be used to store the parameters of the
 * recurrence. In the implementation, it is stored in an object of the class
 * @ref umontreal.ssj.hups.F2wStructure. See the description of this class
 * for more details on how the polynomial is stored.
 *
 * Let @f$\mathbf{x} = (x^{(0)}, …, x^{(p-1)}) \in\mathbb F_2^p@f$ be a
 * @f$p@f$-bit vector. Let us define the function @f$\phi(\mathbf{x}) =
 * \sum_{i=1}^p 2^{-i} x^{(i-1)}@f$. The point set in @f$t@f$ dimensions
 * produced by this class is
 * @f[
 *   \{ (\phi(\mathbf{y}_0),\phi(\mathbf{y}_1),…,\phi(\mathbf{y}_{t-1}): (\mathbf{q}_{0,1},…,\mathbf{q}_{0,r-1})\in\mathbb F_2^{rw}\}
 * @f]
 * where @f$\mathbf{y}_n = (\mathbf{q}_{n,1},…,\mathbf{q}_{n,r})@f$,
 * @f$\mathbf{q}_{n,i}@f$ is the representation of @f$q_{n,i}@f$ under the
 * polynomial basis of @f$\mathbb F_{2^w}@f$ over @f$\mathbb F_2@f$.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class F2wCycleBasedPolyLCG extends CycleBasedPointSetBase2 {

   private F2wStructure param;

   /**
    * Constructs a point set with @f$2^{rw}@f$ points. See the description
    * of the class  @ref umontreal.ssj.hups.F2wStructure for the meaning
    * of the parameters.
    */
   public F2wCycleBasedPolyLCG (int w, int r, int modQ, int step, int nbcoeff,
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
      numBits = param.numBits;
      normFactor = param.normFactor;
      EpsilonHalf = param.EpsilonHalf;
      fillCyclesPolyLCG ();
   }

   /**
    * Constructs a point set after reading its parameters from file
    * `filename`; the parameters are located at line numbered `no` of
    * `filename`. The available files are listed in the description of
    * class  @ref umontreal.ssj.hups.F2wStructure.
    */
   public F2wCycleBasedPolyLCG (String filename, int no)
   {
      param = new F2wStructure (filename, no);
      numBits = param.numBits;
      normFactor = param.normFactor;
      fillCyclesPolyLCG ();
   }


   public String toString ()
   {
       String s = "F2wCycleBasedPolyLCG:" + PrintfFormat.NEWLINE;
       return s + param.toString ();
   }

   private void fillCyclesPolyLCG ()
   {
      int n = 1 << param.getLog2N();
      int i;
      int mask1 = (1 << (31 - param.r * param.w)) - 1;
      int mask2 = ~mask1;
      RandomStream random = new MRG32k3a();
      IntArrayList c;             // Array used to store the current cycle.
      int currentState;           // The state currently visited.

      boolean stateVisited[] = new boolean[n];
      // Indicates which states have been visited so far.
      for (i = 0; i < n; i++)
         stateVisited[i] = false;
      int startState = 0;  // First state of the cycle currently considered.
      numPoints = 0;
      while (startState < n) {
         stateVisited[startState] = true;
         c = new IntArrayList ();
         param.state = startState << param.S;
         c.add (param.state);
         // c.add ((state & mask2) | (mask1 &
         // (random.nextInt(0,2147483647))));
         param.output = param.F2wPolyLCG ();
         // Warning: watch for overflow !!!
         while (param.state != (startState << param.S)) {
            stateVisited[param.state >> param.S] = true;
            // c.add ((param.state&mask2) | (mask1 &
            // (random.nextInt(0,2147483647))));
            c.add (param.state);
            param.output = param.F2wPolyLCG ();
         }
         addCycle (c);
         for (i = startState + 1; i < n; i++)
            if (stateVisited[i] == false)
               break;
         startState = i;
      }
   }
}