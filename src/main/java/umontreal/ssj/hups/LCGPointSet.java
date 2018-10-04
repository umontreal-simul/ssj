/*
 * Class:        LCGPointSet
 * Description:  point set defined via a linear congruential recurrence
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
import cern.colt.list.*;

/**
 * Implements a recurrence-based point set defined via a linear congruential
 * recurrence of the form @f$x_i = a x_{i-1} \mod n@f$ and @f$u_i = x_i /
 * n@f$. The implementation is done by storing the values of @f$u_i@f$ over
 * the set of all cycles of the recurrence.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class LCGPointSet extends CycleBasedPointSet {

      private int a;                      // Multiplier.

   /**
    * Constructs and stores the set of cycles for an LCG with modulus
    * @f$n@f$ and multiplier @f$a@f$. If the LCG has full period length
    * @f$n-1@f$, there are two cycles, the first one containing only 0 and
    * the second one of period length @f$n-1@f$.
    *  @param n            required number of points and modulus of the
    *                      LCG
    *  @param a            generator @f$a@f$ of the LCG
    */
   public LCGPointSet (int n, int a) {
      this.a = a;
      double invn = 1.0 / (double)n;   // 1/n
      DoubleArrayList c;  // Array used to store the current cycle.
      long currentState;  // The state currently visited. 
      int i;
      boolean stateVisited[] = new boolean[n];  
         // Indicates which states have been visited so far.
      for (i = 0; i < n; i++)
         stateVisited[i] = false;
      int startState = 0;    // First state of the cycle currently considered.
      numPoints = 0;
      while (startState < n) {
         stateVisited[startState] = true;
         c = new DoubleArrayList();
         c.add (startState * invn);
         // We use the fact that a "long" has 64 bits in Java.
         currentState = (startState * (long)a) % (long)n;
         while (currentState != startState) {
            stateVisited[(int)currentState] = true;
            c.add (currentState * invn);
            currentState = (currentState * (long)a) % (long)n;
            }
         addCycle (c);
         for (i = startState+1; i < n; i++)
            if (stateVisited[i] == false)
                break;
         startState = i;
         }
      }

   /**
    * Constructs and stores the set of cycles for an LCG with modulus @f$n
    * = b^e + c@f$ and multiplier @f$a@f$.
    */
   public LCGPointSet (int b, int e, int c, int a) {
      this (computeModulus (b, e, c), a);
   }

   private static int computeModulus (int b, int e, int c) {
      int n;
      int i;
      if (b == 2) 
         n = (1 << e);
      else {
         for (i = 1, n = b;  i < e;  i++)  n *= b;
         }
      n += c;
      return n;
   }

   public String toString() {
      StringBuffer sb = new StringBuffer ("LCGPointSet with multiplier a = ");
      sb.append (a);
      return sb.toString();
   }

/**
 * Returns the value of the multiplier @f$a@f$.
 */
public int geta () {
      return a;
   }
}