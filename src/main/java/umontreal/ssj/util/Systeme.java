/*
 * Class:        Systeme
 * Description:  Provides tools related to the system or the computer
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since        January 2011

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
   import java.lang.management.*;
   import java.util.*;
   import java.text.*;
   import java.net.*;

/**
 * This class provides a few tools related to the system or the computer.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class Systeme {
   private Systeme() {}

   /**
    * Returns the name of the host computer.
    *  @return the name of the host computer
    */
   public static String getHostName() {
      String host;
      try {
         InetAddress machine = InetAddress.getLocalHost ();
         host = machine.getHostName ();
      } catch (UnknownHostException uhe) {
         host = "unknown host machine";
      }
      // host = System.getenv("HOSTNAME");
      int j = host.indexOf ('.');
      String name;
      if (j >= 0)
         name = host.substring (0, j);
      else
         name = host;
      return name;
   }

   /**
    * Returns information about the running process: name, id, host name,
    * date and time.
    *  @return information about the running process
    */
   public static String getProcessInfo() {
      StackTraceElement[] stack = Thread.currentThread().getStackTrace ();
      StackTraceElement mai = stack[stack.length - 1];
      String str = mai.getClassName ();

      RuntimeMXBean runtime = ManagementFactory.getRuntimeMXBean();
      Date startTime = new Date(runtime.getStartTime());
      SimpleDateFormat dateFormat =
           new SimpleDateFormat("yyyy-MM-dd' 'HH:mm:ss");

      str += " [" + runtime.getName() + "]";
      str += " [" + dateFormat.format(startTime) + "]";

      return str;
   }

}