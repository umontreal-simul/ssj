/*
 * Class:        Systeme
 * Description:  Provides tools related to the system or the computer
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since        January 2011
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