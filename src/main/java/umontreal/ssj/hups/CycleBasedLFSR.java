/*
 * Class:        CycleBasedLFSR
 * Description:
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
import java.io.*;

/**
 * Linear feedback shift register (LFSR) random number
 * generators @cite rLEC96a, @cite vLEC99a, @cite rNIE92b,
 * produce numbers by generating a sequence of bits from a
 * linear recurrence modulo 2, and forming fractional numbers by taking
 * blocks of successive bits. More precisely, let @f$\mathbb F_2@f$ denote
 * the finite field with two elements (say, 0 and 1). Let @f$P(z) = z^k - a_1
 * z^{k-1} - \cdots- a_k@f$ be a polynomial with coefficients in @f$\mathbb
 * F_2@f$, and consider the recurrence
 * @anchor REF_hups_CycleBasedLFSR_mrg
 * @f[
 *   x_n = a_1 x_{n-1} + \cdots+ a_k x_{n-k}, \tag{mrg}
 * @f]
 * whose characteristic polynomial is @f$P(z)@f$. It should be understood
 * that in ({@link REF_hups_CycleBasedLFSR_mrg mrg}) all
 * computations are performed in @f$\mathbb F_2@f$ (this can be identified
 * with working in integer arithmetic modulo 2). Suppose that
 * @f$\mathbf{s}_0 = (x_0,…,x_{k-1})\in\{0,1\}^k@f$ is fixed and define
 * @anchor REF_hups_CycleBasedLFSR_taus
 * @f[
 *   u_n = \sum_{i=1}^L x_{ns+i-1} 2^{-i}, \tag{taus}
 * @f]
 * where @f$s@f$ and @f$L@f$ are positive integers. If @f$P@f$ is primitive,
 * @f$\mathbf{s}_0\not0@f$, and @f$\rho= 2^k-1@f$ is coprime to @f$s@f$,
 * then the sequences ({@link REF_hups_CycleBasedLFSR_mrg mrg}) 
 * and ({@link REF_hups_CycleBasedLFSR_taus taus}) are both
 * purely periodic with period @f$\rho@f$. Computing @f$u_n@f$ from
 * @f$u_{n-1}@f$ involves performing @f$s@f$ steps of the 
 * recurrence ({@link REF_hups_CycleBasedLFSR_mrg mrg}).
 *
 * Suppose now that we have @f$J@f$ LFSR recurrences, the @f$j@f$-th one
 * having a primitive characteristic polynomial @f$P_j(z)@f$ of degree
 * @f$k_j@f$, and step size @f$s_j@f$. Let @f$\{x_{j,n},  n\ge0\}@f$ be the
 * @f$j@f$-th LFSR sequence,
 * and define @f$x_n = (x_{1,n} + \cdots+ x_{J,n}) \bmod 2@f$ 
 * and @f$u_n@f$ as in ({@link REF_hups_CycleBasedLFSR_taus taus}). Equivalently,
 * if @f$\{u_{j,n},  n\ge0\}@f$ is the output sequence from the @f$j@f$-th
 * LFSR, then @f$u_n = u_{1,n}\oplus\cdots\oplus u_{J,n}@f$ where
 * @f$\oplus@f$ denotes the bitwise exclusive-or in the binary expansion.
 * The sequence @f$\{x_n\}@f$ is called the combined LFSR sequence and a
 * generator that produces this @f$\{u_n\}@f$ is called a *combined LFSR*
 * generator.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class CycleBasedLFSR extends CycleBasedPointSetBase2 {
   private int J = 1; // Nombre de polynomes a combiner
   private int k1 = 0; // Degre du premier polynome.
   private int k2 = 0; // Degre du second polynome.
   private int step1;  // Nombre de repetition de la recurence.
   private int step2;
   private int nbcoeff1; // Nombre de coefficient du polynome.
   private int nbcoeff2;
   private int [] nocoeff1; // Liste des coefficients non nulls.
   private int [] nocoeff2;
   private int [] posa1; // Position des bits de x1 des coefficients a1.
   private int [] posa2;
   private int [] shifta1; // 1 << posa1
   private int [] shifta2;
   private int x1 = 0; // Memorise les x[n-i] de la recurence.
   private int x2 = 0;
   private int state = 0; // Etat present du generateur.
   private int value = 0; // = (x1 ^ x2) & maskX
   private int maskX;  // Recupere les bits 0 a k-1 de x1 ou x2. k = k1+k2
   private int maskX1; // Recupere les bits 0 a k1-1 de x1
   private int maskX2; // Recupere les bits 0 a k2-1 de x2
   private int maskX1Shift; // Recupere les bits k2 a k-1 de x1

   /**
    * Constructs a point set based on a polynomial in base 2 with
    * @f$2^{k_1}@f$ points. `step1` is the number of times the recursion
    * is executed in order to generate the next state. `nbcoeff1` is the
    * number of non-zero coefficients of the polynomial, while `nocoeff1`
    * is the list of the degree of the non-zero coefficients. Example:
    * @f$P(z) = z^9 - z^4 - 1@f$, @f$k_1@f$ = 9, `nbcoeff1` = 3,
    * `nocoeff1` = \{9, 4, 0\}.
    */
   public CycleBasedLFSR (int step1, int nbcoeff1, int [] nocoeff1) {
      this.step1 = step1;
      this.nbcoeff1 = nbcoeff1 - 1;
      this.nocoeff1 = nocoeff1;
      J = 1;
      init ();
   }

   /**
    * Constructs a point set based on a combination of two polynomials in
    * base 2 with @f$2^{k_1 + k_2}@f$ points. The meaning of the
    * parameters is the same as in the case of one polynomial.
    */
   public CycleBasedLFSR (int step1, int step2, int nbcoeff1, int nbcoeff2,
         int [] nocoeff1, int [] nocoeff2) {
      this.step1 = step1;
      this.step2 = step2;
      this.nbcoeff1 = nbcoeff1 - 1;
      this.nbcoeff2 = nbcoeff2 - 1;
      this.nocoeff1 = nocoeff1;
      this.nocoeff2 = nocoeff2;
      J = 2;
      init ();
   }

   /**
    * Constructs a point set after reading its parameters from file
    * `filename`; the parameters associated with number `no` of `filename`
    * corresponds to the <tt>no</tt>-th polynomial. The existing files and
    * the number of polynomials they contain are in the table below. The
    * name of the files describe the number of polynomials @f$J@f$ in the
    * combined LFSR and the number of points @f$2^k@f$ generated. For
    * example, the parameters in file <tt>j1_k11.dat</tt> are based on
    * @f$J = 1@f$ polynomial and generates @f$2^k = 2^{11}@f$ points,
    * while those in file <tt>j2_k17.dat</tt> are based on a combination
    * of @f$J=2@f$ polynomials and generates @f$2^k = 2^{17}@f$ points.
    * Thus to use the 3<em>-th</em> combined LFSR of file
    * <tt>j2_k17.dat</tt>, one must use <tt>CycleBasedLFSR("j2_k17",
    * 3)</tt>.
    *
    * <center>
    *
    * <table class="SSJ-table SSJ-has-hlines">
    * <tr class="bt">
    *   <td class="c bl br">Filename</td>
    *   <td class="c bl br">Num. of polynomials</td>
    * </tr><tr class="bt">
    *   <td class="c bl br"><tt>j1_k11.dat</tt></td>
    *   <td class="c bl br">1</td>
    * </tr><tr>
    *   <td class="c bl br"><tt>j2_k17.dat</tt></td>
    *   <td class="c bl br">6</td>
    * </tr><tr>
    *   <td class="c bl br"><tt>j2_k19.dat</tt></td>
    *   <td class="c bl br">4</td>
    * </tr>
    * </table>
    *
    * </center>
    */
   public CycleBasedLFSR (String filename, int no)
   {
      readFile(filename, no);
      init ();
   }


   // Initialise les variables, et appelle fillCyclesLFSR()
   private void init()
   {
      k1 = nocoeff1[0];
      posa1 = new int[nbcoeff1];
      shifta1 = new int[nbcoeff1];
      for (int i = 0; i < nbcoeff1; ++i) {
         posa1[i] = k1 - nocoeff1[i + 1] - 1;
         shifta1[i] = 1 << posa1[i];
      }

      if (J == 2) {
         k2 = nocoeff2[0];
         posa2 = new int[nbcoeff2];
         shifta2 = new int[nbcoeff2];
         for (int i = 0; i < nbcoeff2; ++i) {
            posa2[i] = k2 - nocoeff2[i + 1] - 1;
            shifta2[i] = 1 << posa2[i];
         }
      }

      if (k1 + k2 > 31) {
         System.err.println("The degree of the combined polynomials must be < 31");
         System.err.println("k1 = " + k1 + " k2 = " + k2);
         System.exit(1);
      }

      maskX1 = (1 << k1) - 1;
      maskX2 = (1 << k2) - 1;
      maskX1Shift = maskX1 << k2;
      maskX = (1 << (k1 + k2)) - 1;

      numBits = k1 + k2;
      normFactor = 1.0 / (1L << numBits);

      fillCyclesLFSR();
   }

   /**
    * This method returns a string containing the polynomials
    * and the stepping parameters.
    */
   public String toString ()
   {
      int i;
      String s = "CycleBasedLFSR:" + PrintfFormat.NEWLINE +
                 "First  Polynome:  Step: " + step1 + "  Coefficients: ";
      for (i = 0; i < nbcoeff1; i++) {
         s += nocoeff1[i] + ", ";
      }
      s += nocoeff1[i];

      if (J == 1) return s;

      s += PrintfFormat.NEWLINE + "Second Polynome:  Step: " +
           step2 + "  Coefficients: ";
      for (i = 0; i < nbcoeff2; i++) {
         s += nocoeff2[i] + ", ";
      }
      s += nocoeff2[i];
      return s;
   }

   // Produit l'etat suivant du generateur
   // x1[n] = a1[1]*x1[n-1] + a1[2]*x1[n-2] + ... + a1[k1]*x1[n-k1]
   // repete step1 fois
   // x2[n] = a2[1]*x2[n-1] + a2[2]*x2[n-2] + ... + a2[k2]*x2[n-k2]
   // repete step2 fois
   // x1[n-1], x1[n-2],...,x1[n-k1] sont les bits 0, 1,..., k1-1, de l'entier x1
   // La position des a1 non nuls est dans le tableau posa1 de taille nbcoeff1
   // shifta1[i] = 1 << posa1[i]
   // value = (x1 ^ x2) & maskX
   // state = x1x2
   private void nextState()
   {
      int i;
      int b = 0;
      int s = 0;
      // Recurence pour le premier polynome.
      for (s = 0; s < step1; s++) {
         b = 0;
         for (i = 0; i < nbcoeff1; i++) {
            b ^= (x1 & shifta1[i]) >> posa1[i];
         }
         x1 = x1 << 1;
         x1 |= b;
      }

      if (J == 1) {
         state = x1 & maskX1;
         value = state;
         return;
      }

      // Recurence pour le second polynome.
      for (s = 0; s < step2; s++) {
         b = 0;
         for (i = 0; i < nbcoeff2; i++) {
            b ^= (x2 & shifta2[i]) >> posa2[i];
         }
         x2 = x2 << 1;
         x2 |= b;
      }

      // value est ce qui est memorise dans les cycles.
      value = (x1 ^ x2) & maskX;

      // L'etat du generateur est constitue de l'etat des recurences x1 et x2.
      // L'etat de la recurence x1 est un nombre de k1 bits situe dans les
      // bits k2 a k-1. k = k1 + k2.
      // L'etat de la recurence x2 est un nombre de k2 bits situe dans les
      // bits k1 a k-1.
      // On forme a partir de ses deux etats un nombre de k bits.
      // Les k1 bits de x1 a gauche, et les k2 bits de x2 a droite.
      // state est utilise pour trouver les cycles.
      // On ne peut pas utiliser value pour trouver les cycles, puisque qu'on
      // ne peut pas retrouver x1 et x2 a partir de x1 ^ x2.
      state = (x1 & maskX1Shift) | ((x2 >>> k1) & maskX2);
   }

   // A partir de l'etat x1x2 on remplit les nombres x1 et x2.
   private void validateState(int x1x2)
   {
      int b = 0;
      int s = 0;
      int i;

      // L'etat de x1 est dans la partie gauche de x1x2
      x1 = x1x2 >> k2;
      // On ramene les k1 bits a la position k2 a k-1
      // Les bits de k a 31 de x1, ne sont jamais utilises.
      for (s = 0; s < k2; s++) {
         b = 0;
         for (i = 0; i < nbcoeff1; i++) {
            b ^= (x1 & shifta1[i]) >> posa1[i];
         }
         x1 = x1 << 1;
         x1 |= b;
      }

      if (J == 1) {
         state = x1 & maskX1;
         value = state;
         return;
      }

      x2 = x1x2 & maskX2;
      // On ramene les k2 bits a la position k1 a k-1
      for (s = 0; s < k1; s++) {
         b = 0;
         for (i = 0; i < nbcoeff2; i++) {
            b ^= (x2 & shifta2[i]) >> posa2[i];
         }
         x2 = x2 << 1;
         x2 |= b;
      }
      // value est ce qui est memorise dans les cycles.
      value = (x1 ^ x2) & maskX;
   }

   // Remplit les cycles avec le generateur LFSR.
   // Le nombre de valeurs dans l'ensemble des cycles est de 2^(k1+k2).
   // Cette methode fonctionne pour tous les polynomes, meme ceux dont
   // la periode n'est pas maximale.
   private void fillCyclesLFSR ()
   {
      int n = 1 << (k1 + k2);  // Nombre de points dans l'ensemble.
      IntArrayList c;          // Array used to store the current cycle.
      int i;
      boolean stateVisited[] = new boolean[n];

      // Indicates which states have been visited so far.
      for (i = 0; i < n; i++)
         stateVisited[i] = false;
      int startState = 0;  // First state of the cycle currently considered.
      numPoints = 0;
      while (startState < n) {
         stateVisited[startState] = true;
         c = new IntArrayList ();
         c.add (value);
         nextState();
         while (state != startState) {
            stateVisited[state] = true;
            c.add (value);
            nextState();
         }
         addCycle (c);
//         System.out.println("Size of Cycle: " + c.size());
         for (i = startState + 1; i < n; i++)
            if (stateVisited[i] == false)
               break;
         startState = i;
         validateState(i);
      }
   }

   // Lit le fichier contenant les polynomes pour plusieurs generateurs.
   // no indique le choix du generateur dans le fichier.
   private void readFile(String fileName, int no) {
      BufferedReader input = null;
      try {
         if ((new File (fileName)).exists()) {
            input = new BufferedReader (new FileReader (fileName));
         } else {
            DataInputStream dataInput;
            dataInput = new DataInputStream (
               F2wStructure.class.getClassLoader().getResourceAsStream (
                   "umontreal/ssj/hups/dataLFSR/" + fileName));
            input = new BufferedReader (new InputStreamReader (dataInput));
         }
      }
      catch (FileNotFoundException e) {
         System.err.println ("File " + fileName + " not found" +
                             PrintfFormat.NEWLINE );
         System.exit(1);
      }

      // La numerotation des generateurs debute a 1
      if (no < 1) no = 1;

      int [] numbers;

      // La premiere ligne valide est le nombre de polynomes
      String line = readOneLine(input);
      numbers = lineToNumbers(line);
      J = numbers[0];
      if (J != 1  &&  J != 2) {
         System.err.println("Error: J = " + J + PrintfFormat.NEWLINE +
           "CycleBasedLFSR works only for the cases of one or two polynomials");
         System.exit(1);
      }

      // Nombre ne lignes qu'il faut sauter avant de lire les polynomes.
      int nbLines = (J + 1) * (no - 1);

      // On saute nbLines valides
      for (int i = 0; i < nbLines; i++) {
         line = readOneLine(input);
         if (line == null) {
            // Il n'y a aucune ligne qui correspond a no.
            System.err.println("Error CycleBasedLFSR:" +
                   PrintfFormat.NEWLINE + " no data in file " +
                  fileName + " for " + no + "-th LFSR");
            System.exit(1);
         }
      }

      // Lecture des valeurs des J steps.
      line = readOneLine(input);
      if (line == null) {
         // Il n'y a aucune ligne qui correspond a no.
          System.err.println("Error CycleBasedLFSR:" +
                  PrintfFormat.NEWLINE + " no data in file " +
                  fileName + " for " + no + "-th LFSR");
         System.exit(1);
      }
      numbers = lineToNumbers(line);

      // Lecture des J polynomes.
      if (J == 1) {
         step1 = numbers[0];
         line = readOneLine(input);
         nocoeff1 = lineToNumbers(line);
         nbcoeff1 = nocoeff1.length - 1;
      }
      else if (J == 2) {
         step1 = numbers[0];
         step2 = numbers[1];
         line = readOneLine(input);
         nocoeff1 = lineToNumbers(line);
         nbcoeff1 = nocoeff1.length - 1;
         line = readOneLine(input);
         nocoeff2 = lineToNumbers(line);
         nbcoeff2 = nocoeff2.length - 1;
      }
   }

   // Lit une ligne du fichier en ignorant les lignes blanches
   // et les commentaires.
   private String readOneLine(BufferedReader input)
   {
      String line;

      try {
         while (true) {
            line = input.readLine();
            if (line == null) return null;
            line = line.trim();

            // On saute les lignes blanches.
            if (line.length() == 0) {
               line = input.readLine();
               continue;
            }

            // On saute les lignes de commentaires
            if (line.charAt(0) == '#') {
               line = input.readLine();
               continue;
            }

            // On elimine les commentaires a la fin de la ligne
            int index = line.indexOf('#');
            if (index >= 0) {
               line = line.substring(0, index).trim();
            }
            return line;
         }
      }
      catch (IOException e) {
         System.err.println(e);
         System.exit(1);
      }
      return null;
   }

   // Convertit la String line en tableau d'entiers
   private int [] lineToNumbers(String line)
   {
      int [] numbers;
      String [] snumbers;

      // On separent les nombres.
      snumbers = line.split("\\s++");

      int nb = snumbers.length;
      numbers = new int[nb];
      for (int i = 0; i < nb; i++) {
         numbers[i] = Integer.valueOf(snumbers[i]).intValue();
      }
      return numbers;
   }

}