/*
 * Class:        ThreadCPUTimeChrono
 * Description:  Compute the CPU time for a single thread
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

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;

/**
 * Extends the  @ref AbstractChrono class to compute the CPU time for a
 * single thread. It is available only under Java 1.5 which provides
 * platform-independent facilities to get the CPU time for a single thread
 * through management API.
 *
 * Note that this chrono might not work properly on some systems running
 * Linux because of a bug in Sun’s implementation or Linux kernel. For
 * instance, this class unexpectedly computes the global CPU time under
 * Fedora Core&nbsp;4, kernel 2.6.17 and JRE version 1.5.0-09. With Fedora
 * Core&nbsp;6, kernel 2.6.20, the function is working properly. As a result,
 * one should not rely on this bug to get the global CPU time.
 *
 * Note that the above bug does not prevent one from using this chrono to
 * compute the CPU time for a single-threaded application. In that case, the
 * global CPU time corresponds to the CPU time of the current thread.
 *
 * Running timer fonctions when the associated thread is dead will return 0.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class ThreadCPUTimeChrono extends AbstractChrono {
   private  long           myThreadId;
   static   ThreadMXBean   threadMXBean = null;

   protected void getTime (long[] tab) {
      long rawTime = getTime();
      final long DIV = 1000000000L;
      long seconds = rawTime/DIV;
      long micros = (rawTime % DIV)/1000L;
      tab[0] = seconds;
      tab[1] = micros;
   }

   protected long getTime() {
      if (threadMXBean == null) {
         // We use lazy initialization to avoid a potential exception being wrapped into a confusing
         // ExceptionInInitializerError. That would happen if this initialization was in a static block instead of
         // in this method.
         threadMXBean = ManagementFactory.getThreadMXBean();
         if (!threadMXBean.isThreadCpuTimeEnabled())
            // Call this only when necessary, because this can throw a SecurityException if
            // run under a security manager.
            threadMXBean.setThreadCpuTimeEnabled (true);
      }
      long time = threadMXBean.getThreadCpuTime(myThreadId);
      return time < 0 ? 0 : time;
   }

   /**
    * Constructs a `ThreadCPUTimeChrono` object associated with current
    * thread and initializes it to zero.
    */
   public ThreadCPUTimeChrono() {
      super();
      myThreadId = Thread.currentThread().getId();
      init();
   }

   /**
    * Constructs a `ThreadCPUTimeChrono` object associated with the given
    * Thread variable and initializes it to zero.
    */
   public ThreadCPUTimeChrono(Thread inThread) {
      super();
      myThreadId = inThread.getId();
      init();
   }

}