/*
 * Class:        TimeUnit
 * Description:  Represents a time unit for conversion of time durations.
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
package umontreal.ssj.util;

/**
 * Represents a time unit for conversion of time durations. A time unit
 * instance can be used to get information about the time unit and as a
 * selector to perform conversions. Each time unit has a short name used when
 * representing a time unit, a full descriptive name, and the number of hours
 * corresponding to one unit.
 *
 * <div class="SSJ-bigskip"></div>
 */
public enum TimeUnit {

   /**
    * @name enum values
    * @{
    */

   /**
    * Represents a nanosecond which has short name `ns`.
    */
   NANOSECOND

("ns", "nanosecond", 1.0/3600000000000.0),

   /**
    * Represents a microsecond which has short name `us`.
    */
   MICROSECOND

("us", "microsecond", 1.0/3600000000.0),

   /**
    * Represents a millisecond which has short name `ms`.
    */
   MILLISECOND

("ms", "millisecond", 1.0/3600000.0),

   /**
    * Represents a second which has short name `s`.
    */
   SECOND

("s", "second", 1.0/3600.0),

   /**
    * Represents a minute which has short name `min`.
    */
   MINUTE

("min", "minute", 1.0/60.0),

   /**
    * Represents an hour which has short name `h`.
    */
   HOUR

("h", "hour", 1.0),

   /**
    * Represents a day which has short name `d`.
    */
   DAY

("d", "day", 24.0),

   /**
    * Represents a week which has short name `w`.
    */
   WEEK

("w", "week", 24.0*7);


   private String shortName;
   private String longName;
   private transient double numHours;

   private TimeUnit (String shortName, String longName, double numHours) {
      this.shortName = shortName;
      this.longName = longName;
      this.numHours = numHours;
   }

   /**
    * @}
    */

   /**
    * Returns the short name representing this unit in a string specifying
    * a time duration.
    *  @return the short name of this time unit.
    */
   public String getShortName() {
      return shortName;
   }

   /**
    * Returns the long name of this time unit.
    *  @return the long name of this time unit.
    */
   public String getLongName() {
      return longName;
   }

   /**
    * Calls  #getLongName.
    *  @return the result of  #getLongName.
    */
   public String toString() {
      return longName;
   }

   /**
    * Returns this time unit represented in hours. This returns the number
    * of hours corresponding to one unit.
    *  @return the time unit represented in hours.
    */
   public double getHours() {
      return numHours;
   }

   /**
    * Converts `value` expressed in time unit `srcUnit` to a time duration
    * expressed in `dstUnit` and returns the result of the conversion.
    *  @param value        the value being converted.
    *  @param srcUnit      the source time unit.
    *  @param dstUnit      the destination time unit.
    *  @return the converted value.
    *
    *  @exception NullPointerException if `srcUnit` or `dstUnit` are
    * `null`.
    */
   public static double convert (double value, TimeUnit srcUnit,
                                 TimeUnit dstUnit) {
      double hours = value*srcUnit.getHours();
      return hours/dstUnit.getHours();
   }
}