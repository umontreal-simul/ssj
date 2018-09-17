package umontreal.ssj.rng;

import umontreal.ssj.util.PrintfFormat;

/*
 * Class:        F2wPoly
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author:
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


class F2wPoly {
   private int r;                  //dimension
   private F2w.F2wElem coeff[];    //coefficients non-nuls
   private int noCoeff[];          //position des coefficients non-nuls
   private F2w f2wBase;            //l'espace des F2w

   private F2wPolyElem z_i[];      //resultats du calculs des z^i dans F2wPoly


   /*
     Constructeur :
     Les coefficients coeff sont ceux du polynome P(z).
     Ils doivent avoir ete cree par la base f2wBase.
   */
   public F2wPoly(int r, F2w.F2wElem coeff[], int noCoeff[], F2w f2wBase) {
      this.r = r;

      this.coeff = new F2w.F2wElem[coeff.length];
      this.noCoeff = new int[noCoeff.length];
      for(int i = 0; i < coeff.length; i++) {
         this.coeff[i] = coeff[i];
         this.noCoeff[i] = noCoeff[i];
      }

      this.f2wBase = f2wBase;


      //initialisation des z (pre-calculs)
      z_i = new F2wPolyElem[2 * r];
      z_i[0] = new F2wPolyElem();
      for(int i = 0; i < r; i++)
         z_i[0].value[i].value = 0;
      z_i[0].value[0].value = 0x80000000;

      for(int i = 1; i < 2*r; i++)
         z_i[i] = z_i[i-1].multiplyZ();
   }


   /*
     Differentes methodes permettant d'acceder aux constructeurs de
     F2wPolyElem (on ne peut y acceder directement puisque F2wPolyElem
     a besoin de connaitre a quel F2wPoly il appartient).
   */
   public F2wPolyElem createElem() {
      return new F2wPolyElem();
   }
   public F2wPolyElem createElem(F2wPolyElem that) {
      return new F2wPolyElem(that);
   }
   public F2wPolyElem createElem(int[] value) {
      return new F2wPolyElem(value);
   }

   /*
     Retourne le F2wPolyElem qui est egal a au polynome "z"
   */
   public F2wPolyElem createZ() {
      int[] val = new int[r];
      for(int i = 0; i < r; i++)
         val[i] = 0;
      val[1] = 0x80000000;
      return createElem(val);
   }



   /*
     Represente l'espace (F2w) dont font partie les coefficients des
     polynomes.

     F2w est definie comme etant F2[zeta]\Q(zeta).

     Note : la classe a ete concue pour que w = 32.
   */
   public static class F2w {
      private int w;       //dimension
      private int modQ;    //modulo

      private F2wElem[] zeta_i;


      public F2w(int modQ) {
         this(32, modQ);
      }

      private F2w(int w, int modQ) {
         this.w = w;
         this.modQ = modQ;

         // initialisation des zetas
         zeta_i = new F2wElem[w];
         zeta_i[0] = new F2wElem(0x80000000);
         for(int i = 1; i < w; i++)
            zeta_i[i] = zeta_i[i-1].multiplyZeta(2);
      }


      /*
        Methodes permettant d'acceder aux constructeurs de F2wElem.


      */
      public F2wElem createElem() {
         return new F2wElem();
      }
      public F2wElem createElem(F2wElem that) {
         return new F2wElem(that);
      }
      public F2wElem createElem(int val) {
         return new F2wElem(val);
      }

      public int getDim() {
         return w;
      }
      public int getModulo() {
         return modQ;
      }


      /*
        Classe representant les elements de F2w.
        Le bit le plus significatif (celui en 0x80000000) represente
        zeta^0 et le moins significatif (celui en 0x00000001) represente
        zeta^31. La multiplication par zeta se fait donc par un rigth-shift,
        suivis d'un modulo si le resultat depasse zeta^32.
      */
      public class F2wElem {
         private int value;

         //constructeurs
         private F2wElem() {
            value = 0;
         }
         private F2wElem(F2wElem that) {
            if(this.getBase() != that.getBase())
               throw new IllegalArgumentException
               ("The copied F2wElem must originate from the same F2w.");
            this.value = that.value;
         }
         private F2wElem(int val) {
            value = val;
         }

         public F2w getBase() {
            return F2w.this;
         }
         public int getValue() {
            return value;
         }

         // calcule this * zeta^k dans F2w
         public F2wElem multiplyZeta(int k) {
            int res = value;

            if(k == 0)
               return new F2wElem(res);
            else {
               for(int i = 0; i < k; i++)
                  if((1 & res) != 0)
                     res = (res >>> 1) ^ modQ;
                  else
                     res >>>= 1;
               return new F2wElem(res);
            }
         }

         // calcule this * that dans F2w
         public F2wElem multiply(F2wElem that) {
            if(this.getBase() != that.getBase())
               throw new IllegalArgumentException
               ("Both F2wElem must originate from the same F2w.");

            int res = 0;
            int verif = 1;

            for(int i = w - 1; i >= 0; i--) {
               if((that.value & verif) != 0)
                  res ^= this.multiplyZeta(i).value;
               verif <<= 1;
            }

            return new F2wElem(res);
         }

         // calcule (this)^2 dans F2w
         public F2wElem square() {
            int res = 0;

            for(int i = 0; i < w; i++)
               if((value & (0x80000000 >>> i)) != 0)
                  res ^= zeta_i[i].value;

            return new F2wElem(res);
         }

         public String toString() {
            StringBuffer sb = new StringBuffer();

            int temp = value;

            for(int i = 0; i < 32; i++) {
               sb.append((temp & 1) == 1 ? '1' : '0');
               temp >>>= 1;
            }

            return sb.reverse().toString();
         }
      }

   }



   /*
     Represente un polynome faisant partie de F2wPoly.

     Le premier element de value est le coefficient de z^0, tandis que
     le dernier est le coefficient de z^(r-1).
   */
   public class F2wPolyElem {
      private F2w.F2wElem[] value;

      private F2wPolyElem() {
         value = new F2w.F2wElem[r];
         for(int i = 0; i < r; i++)
            value[i] = f2wBase.createElem();
      }

      private F2wPolyElem(F2wPolyElem that) {
         if(this.getBase() != that.getBase())
            throw new IllegalArgumentException
            ("The copied F2wPolyElem must come from the same F2wPoly.");

         value = new F2w.F2wElem[r];
         for(int i = 0; i < r; i++)
            this.value[i] = f2wBase.createElem(that.value[i]);
      }

      private F2wPolyElem(int[] value) {
         if(r != value.length)
            throw new IllegalArgumentException
            ("Array length must be equal to r (" + r + ")");

         this.value = new F2w.F2wElem[r];
         for(int i = 0; i < r; i++)
            this.value[i] = f2wBase.createElem(value[i]);
      }


      public F2wPoly getBase() {
         return F2wPoly.this;
      }


      //multiplie par this par z^1 dans F2w[z]
      public F2wPolyElem multiplyZ() {
         F2wPolyElem res = new F2wPolyElem();

         res.value[0].value = 0;
         for(int i = 1; i < r; i++)
            res.value[i].value = this.value[i-1].value;

         for(int i = 0; i < noCoeff.length; i++)
            res.value[noCoeff[i]].value ^=
               this.value[r-1].multiply(coeff[i]).value;

         return res;
      }

      //calcule this * that dans F2w[z]
      public F2wPolyElem multiply(F2wPolyElem that) {
         if(this.getBase() != that.getBase())
            throw new IllegalArgumentException
            ("Both F2wPolyElem must originate from the same F2wPoly.");

         F2wPolyElem res = new F2wPolyElem();

         for(int i = 0; i < r; i++)
            for(int j = 0; j < r; j++) {
               F2w.F2wElem temp = this.value[i].multiply(that.value[j]);
               for(int k = 0; k < r; k++)
                  res.value[k].value ^=
                     z_i[i+j].value[k].multiply(temp).value;
            }
         return res;
      }


      //calcule this^(2^d)
      public F2wPolyElem exponentiateBase2 (int d) {
         F2wPolyElem res = new F2wPolyElem(this);
         F2wPolyElem temp= new F2wPolyElem();

         for (int i = 0; i < d; i++) {
            for(int j = 0; j < r; j++)
               temp.value[j].value = 0;

            for(int j = 0; j < r; j++) {
               F2w.F2wElem coeff = res.value[j].square();
               for(int k = 0; k < r; k++)
                  temp.value[k].value ^=
                     z_i[2*j].value[k].multiply(coeff).value;
            }

            for(int j = 0; j < r; j++)
               res.value[j].value = temp.value[j].value;
         }

         return res;
      }


      public void copyFrom(F2wPolyElem that) {
         if(this.getBase() != that.getBase())
            throw new IllegalArgumentException
            ("Both F2wPolyElem must originate from the same F2wPoly.");

         for(int i = 0; i < r; i++)
            this.value[i].value = that.value[i].value;
      }
      public void copyFrom(int[] val) {
         if(r != value.length)
            throw new IllegalArgumentException
            ("Array length must be equal to r (" + r + ")");

         for(int i = 0; i < r; i++)
            this.value[i].value = val[i];
      }
      public void copyTo(int[] val) {
         if(r != value.length)
            throw new IllegalArgumentException
            ("Array length must be equal to r (" + r + ")");

         for(int i = 0; i < r; i++)
            val[i] = this.value[i].value;
      }


      public String toString() {
         StringBuffer sb = new StringBuffer("{");

         for(int i = 0; i < r - 1; i++)
            sb.append(value[i].toString() + ", " +
                PrintfFormat.NEWLINE + " ");
         if(r > 0)
            sb.append(value[r-1].toString());
         sb.append("}");

         return sb.toString();
      }

   }


   public static void main(String[] args) {
      F2w f1 = new F2w(32, 0x00010002);
      F2w f2 = new F2w(32, 0x10204080);

      F2w.F2wElem e1 = f1.createElem(0x12345678);
      F2w.F2wElem e2 = f2.createElem(0x19414111);
      F2w.F2wElem e3 = e1.multiply(e2);

      /*
      F2w f2w = new F2w(32, 0xFA4F9B3F);

      F2wPoly poly = new F2wPoly(25,
                                 new F2w.F2wElem[]{f2w.createElem(0xE6A68D20),
                                                   f2w.createElem(0x287AB842)},
                                 new int[]{7,
                                           0},
                                 f2w);

      F2wPolyElem gen = poly.createElem
         (new int[]{0x95F24DAB, 0x0B685215, 0xE76CCAE7, 0xAF3EC239,
                    0x715FAD23, 0x24A590AD, 0x69E4B5EF, 0xBF456141,
                    0x96BC1B7B, 0xA7BDF825, 0xC1DE75B7, 0x8858A9C9,
                    0x2DA87693, 0xB657F9DD, 0xFFDC8A8F, 0x8121DA71,
                    0x8B823ECB, 0x885D05F5, 0x4E20CD47, 0x5A9AD5D9,
                    0x512C0C03, 0xEA857CCD, 0x4CC1D30F, 0x8891A8A1,
                    0xA6B7AADB});


      for(int i = 0; i < poly.z_i.length; i++)
         {
            System.out.println(i + " :");
            System.out.println(poly.z_i[i]);
            System.out.println();
         }
      */

   }

}

