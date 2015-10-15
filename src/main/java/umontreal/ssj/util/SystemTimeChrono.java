/*
 * Class:        SystemTimeChrono
 * Description:  
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Éric Buist
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
package umontreal.ssj.util;

/**
 * Extends the  @ref AbstractChrono class to compute the total system time
 * using Java’s builtin `System.nanoTime`. The system can be used as a rough
 * approximation of the CPU time taken by a program if no other tasks are
 * executed on the host while the program is running.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class SystemTimeChrono extends AbstractChrono {

   protected void getTime (long[] tab) {
      long rawTime = System.nanoTime();
      final long DIV = 1000000000L;
      long seconds = rawTime/DIV;
      long micros = (rawTime % DIV)/1000L;
      tab[0] = seconds;
      tab[1] = micros;
   }

   /**
    * Constructs a new chrono object and initializes it to zero.
    */
   public SystemTimeChrono() {
      super();
      init();
   }

}