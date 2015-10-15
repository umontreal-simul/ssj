/*
 * Class:        LCGPointSet
 * Description:  point set defined via a linear congruential recurrence
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
      StringBuffer sb = new StringBuffer ("LCGPointSet:" +
                                           PrintfFormat.NEWLINE);
      sb.append (super.toString());
      sb.append (PrintfFormat.NEWLINE + "Multiplier a: ");
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