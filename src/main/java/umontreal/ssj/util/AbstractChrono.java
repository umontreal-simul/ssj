/*
 * Class:        AbstractChrono
 * Description:  calculates CPU time of parts of a program
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
package umontreal.ssj.util;

/**
 * `AbstractChrono` is a class that acts as an interface to the system clock
 * and calculates the CPU or system time consumed by parts of a program.
 *
 * Every object of class `AbstractChrono` acts as an independent *stopwatch*.
 * Several `AbstractChrono` objects can run at any given time. The method
 * #init resets the stopwatch to zero,  #getSeconds,  #getMinutes and
 * #getHours return its current reading, and  #format converts this reading
 * to a  String. The returned value includes the execution time of the method
 * from `AbstractChrono`.
 *
 * Below is an example of how it may be used. A stopwatch named `timer` is
 * constructed (and initialized). When 2.1 seconds of CPU time have been
 * consumed, the stopwatch is read and reset to zero. Then, after an
 * additional 330 seconds (or 5.5 minutes) of CPU time, the stopwatch is read
 * again and the value is printed to the output in minutes.
 *
 * @code
 *
 * AbstractChrono timer = new Chrono();
 *
 * @endcode
 *
 * @code
 *
 *     \vdots
 * @endcode
 *   (*suppose 2.1 CPU seconds are used here*.)
 *
 * @code
 *
 * double t = timer.getSeconds();               // Here, t = 2.1
 * timer.init();
 * t = timer.getSeconds();                      // Here, t = 0.0
 *
 * @endcode
 *
 * @code
 *     \vdots
 * @endcode
 *   (*suppose 330 CPU seconds are used here*.)
 *
 * @code
 *
 * t = timer.getMinutes();                      // Here, t = 5.5
 * System.out.println (timer.format());         // Prints: 0:5:30.00
 *
 * @endcode
 *
 * <strong>Warning:</strong><br>
 *
 *  Even though the ANSI/ISO macro <tt>CLOCKS_PER_SEC = 1000000</tt> is the
 * number of clock ticks per second for the value returned by the `clock`
 * function (so this function returns the number of microseconds), on some
 * systems the value returned by `clock` wraps around to 0 after about 36
 * minutes when the type `long` used to measure time has only 32 bits (on
 * some systems, the clock may wrap around to 0 after about 72 minutes if
 * time is measured using the 32-bit type <tt>unsigned long</tt>). When the
 * macro <tt>USE_ANSI_CLOCK</tt> in module `gdef` is undefined, a non-ANSI
 * clock is used. It calls the POSIX function `times` to get the CPU time
 * used by a program and the non-ANSI macro <tt>CLK_TCK</tt> is used to get
 * the real number of clock ticks per second for the system.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public abstract class AbstractChrono {

   private long m_second;
   private long m_microsec;
   private long[] now = new long[2];
   // tab[0] = seconds, tab[1] = microseconds

   protected abstract void getTime (long[] tab);

/**
 * @name Timing functions
 * @{
 */
public AbstractChrono() {
   }

/**
 * Initializes this `AbstractChrono` to zero.
 */
public void init() {
      getTime (now);
      m_second = now[0];
      m_microsec = now[1];
   }

   /**
    * Returns the CPU time in seconds used by the program since the last
    * call to  #init for this `AbstractChrono`.
    *  @return the number of seconds
    */
   public double getSeconds() {
      getTime (now);
      double time = (now[1] - m_microsec)/1000000.0
             + (now[0] - m_second);
      return time;
   }

   /**
    * Returns the CPU time in minutes used by the program since the last
    * call to  #init for this `AbstractChrono`.
    *  @return the number of minutes
    */
   public double getMinutes() {
      getTime (now);
      double time = (now[1] - m_microsec)/1000000.0
             + (now[0] - m_second);
      return time*1.666666667*0.01;
   }

   /**
    * Returns the CPU time in hours used by the program since the last
    * call to  #init for this `AbstractChrono`.
    *  @return the number of hours
    */
   public double getHours() {
      getTime (now);
      double time = (now[1] - m_microsec)/1000000.0
             + (now[0] - m_second);
      return time*2.777777778*0.0001;
   }

   /**
    * Converts the CPU time used by the program since its last call to
    * #init for this <tt>AbstractChrono</tt> to a  String in the
    * `HH:MM:SS.xx` format.
    *  @return the string representation of the CPU time
    */
   public String format() {
      return format (getSeconds());
   }

   /**
    * Converts the time `time`, given in seconds, to a  String in the
    * `HH:MM:SS.xx` format.
    *  @return the string representation of the time `time`
    */
   public static String format (double time) {
      int second, hour, min, centieme;
      hour = (int)(time/3600.0);
      if (hour > 0) time -= ((double)hour*3600.0);
      min = (int)(time/60.0);
      if (min > 0) time -= ((double)min*60.0);
      second = (int)time;
      centieme = (int)(100.0*(time - (double)second) + 0.5);
      return String.valueOf (hour) + ":" +
                      min + ":" +
                      second + "." +
                      centieme;
   }

}

/**
 * @}
 */