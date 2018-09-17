/*
 * Class:        DistributionFactory
 * Description:  allows the creation of distribution objects from a string
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
package umontreal.ssj.probdist;
import java.lang.reflect.*;
import java.util.StringTokenizer;

/**
 * This class implements a string API for the package `probdist`. It uses
 * Java Reflection to allow the creation of probability distribution objects
 * from a string. This permits one to obtain distribution specifications from
 * a file or dynamically from user input during program execution. This
 * string API is similar to that of
 * [UNURAN](http://statistik.wu-wien.ac.at/unuran/)
 *  @cite iLEY02a&thinsp;.
 *
 * The (static) methods of this class invoke the constructor specified in the
 * string. For example,
 *
 * @code
 *
 *  d = DistributionFactory.getContinuousDistribution ("NormalDist (0.0,
 * 2.5)");
 *
 * @endcode
 *
 *  is equivalent to
 *
 * @code
 *
 *  d = NormalDist (0.0, 2.5);
 *
 * @endcode
 *
 * The string that specifies the distribution (i.e., the formal parameter
 * `str` of the methods) must be a valid call of the constructor of a class
 * that extends  @ref ContinuousDistribution or  @ref DiscreteDistribution,
 * and all parameter values must be numerical values (variable names are not
 * allowed).
 *  If no parentheses follow, the default parameters are used, i.e. a
 * no-parameter constructor is searched and invoked. Distribution parameters
 * are surrounded with parentheses and are separated from each other using
 * commas. The order and types of the parameters should be the same as the
 * corresponding class constructor. When specifying decimal values, the dot
 * should be used for decimal separation, not the comma. When specifying a
 * `float` or a `double`, one can use `-infinity` or `infinity` to denote
 * @f$-\infty@f$ and @f$+\infty@f$, respectively. This is parsed to
 * <tt>Double.NEGATIVE_INFINITY</tt> and <tt>Double.POSITIVE_INFINITY</tt>,
 * respectively. However, this is not accepted by all probability
 * distributions.
 *
 * For example, if one uses the string `Normal`, a  @ref NormalDist object
 * with @f$\mu= 0@f$ and @f$\sigma= 1@f$ will be created. If one uses
 * `Exponential(2.5)`, an  @ref ExponentialDist object with
 * @f$\lambda=2.5@f$ will be constructed.
 *
 * The distribution parameters can also be estimated from a set of
 * observations instead of being passed to the constructor. In that case, one
 * passes the vector of observations, and the constructor estimates the
 * parameters by the maximum likelihood method.
 *
 * <div class="SSJ-bigskip"></div>
 *
 * @ingroup probdist_general
 */
public class DistributionFactory {
   private DistributionFactory() {}   //  ????   Utile?

   public static Distribution getDistribution (String str) {
      // Extracts the name of the distribution.
      // If there is an open parenthesis, the name contains all the 
      // non-space characters preceeding it.If not,the name is the full string.

      int i = 0;
      str = str.trim();

      int idx = str.indexOf ('(', i);
      String distName;
      if (idx == -1)
         distName = str.substring (i).trim();
      else
         distName = str.substring (i, idx).trim();
 
      // Try to find the class in probdist package.
      Class<?> distClass;
      if (distName.equals ("String"))
         throw new IllegalArgumentException ("Invalid distribution name: " 
                                             + distName);
      try {
         distClass = Class.forName ("umontreal.ssj.probdist." 
                                    + distName);
      }
      catch (ClassNotFoundException e) {
         // Look for a fully qualified classname whose constructor
         //  matches this string.
         try {
            distClass = Class.forName (distName);
            // We must check if the class implements Distribution 
            if (Distribution.class.isAssignableFrom(distClass) == false)
               throw new IllegalArgumentException 
                  ("The given class is not a Probdist distribution class.");
         }
         catch (ClassNotFoundException ex) {
            throw new IllegalArgumentException ("Invalid distribution name: " 
                                                + distName);
         }
      }

      String paramStr = "";
      if (idx != -1) {
         // Get the parameters from the string.
         int parFrom = idx;
         int parTo = str.lastIndexOf (')');
         // paramStr will contain the parameters without parentheses.
         paramStr = str.substring (parFrom + 1, parTo).trim();
         if (paramStr.indexOf ('(') != -1 || paramStr.indexOf (')') != -1)
            //All params are numerical,so parenthesis nesting is forbidden here
            throw new IllegalArgumentException ("Invalid parameter string: " 
                                                + paramStr);
      }

      if (paramStr.equals ("")) {
         // No parameter is given to the constructor.
         try {
            return (Distribution) distClass.newInstance();
         }
         catch (IllegalAccessException e) {
            throw new IllegalArgumentException 
                                         ("Default parameters not available");
         }
         catch (InstantiationException e) {
            throw new IllegalArgumentException 
                                         ("Default parameters not available");
         }
      }

      // Find the number of parameters and try to find a matching constructor.
      // Within probdist, there are no constructors with the same
      // number of arguments but with different types.
      // This simplifies the constructor selection scheme.
      StringTokenizer paramTok = new StringTokenizer (paramStr, ",");
      int nparams = paramTok.countTokens();
      Constructor[] cons = distClass.getConstructors();
      Constructor distCons = null;
      Class[] paramTypes = null;
      // Find a public constructor with the correct number of parameters.
      for (i = 0; i < cons.length; i++) {
         if (Modifier.isPublic (cons[i].getModifiers()) &&
             ((paramTypes = cons[i].getParameterTypes()).length == nparams)) {
            distCons = cons[i];
            break;
         }
      }
      if (distCons == null)
         throw new IllegalArgumentException ("Invalid parameter number");

      // Create the parameters for the selected constructor.
      Object[] instParams = new Object[nparams];
      for (i = 0; i < nparams; i++) {
         String par = paramTok.nextToken().trim();
         try {
            // We only need a limited set of parameter types here.
            if (paramTypes[i] == int.class)
               instParams[i] = new Integer (par);
            else if (paramTypes[i] == long.class)
               instParams[i] = new Long (par);
            else if (paramTypes[i] == float.class) {
               if (par.equalsIgnoreCase ("infinity") || par.equalsIgnoreCase 
                                                                 ("+infinity"))
                  instParams[i] = new Float (Float.POSITIVE_INFINITY);
               else if (par.equalsIgnoreCase ("-infinity"))
                  instParams[i] = new Float (Float.NEGATIVE_INFINITY);
               else
                  instParams[i] = new Float (par);
            }
            else if (paramTypes[i] == double.class) {
               if (par.equalsIgnoreCase ("infinity") || par.equalsIgnoreCase
                                                                 ("+infinity"))
                  instParams[i] = new Double (Double.POSITIVE_INFINITY);
               else if (par.equalsIgnoreCase ("-infinity"))
                  instParams[i] = new Double (Double.NEGATIVE_INFINITY);
               else
                  instParams[i] = new Double (par);
            }
            else
               throw new IllegalArgumentException
                  ("Parameter " + (i+1) + " type " + paramTypes[i].getName() +
                   "not supported");
         }
         catch (NumberFormatException e) {
            throw new IllegalArgumentException
               ("Parameter " + (i+1) + " of type " +
                paramTypes[i].getName()+" could not be converted from String");
         }
      }

      // Try to instantiate the distribution class.
      try {
         return (Distribution) distCons.newInstance (instParams);
      }
      catch (IllegalAccessException e) {
         return null;
      }
      catch (InstantiationException e) {
         return null;
      }
      catch (InvocationTargetException e) {
         return null;
      }
   }

/**
 * Uses the Java Reflection API to construct a  @ref ContinuousDistribution ,
 * @ref DiscreteDistributionInt or  @ref DiscreteDistribution object, as
 * specified by the string `str`. This method throws exceptions if it cannot
 * parse the given string and returns `null` if the distribution object
 * simply could not be created due to a Java-specific instantiation problem.
 *  @param str          distribution specification string
 *  @return a distribution object or `null` if it could not be instantiated
 *
 *  @exception IllegalArgumentException if parsing problems occured when
 * reading `str`
 */
 /**
 * Uses the Java Reflection API to construct a  @ref ContinuousDistribution
 * object by estimating parameters of the distribution using the maximum
 * likelihood method based on the @f$n@f$ observations in table @f$x[i]@f$,
 * @f$i = 0, 1, …, n-1@f$.
 *  @param distName     the name of the distribution to instanciate
 *  @param x            the list of observations to use to evaluate
 *                      parameters
 *  @param n            the number of observations to use to evaluate
 *                      parameters
 */
@SuppressWarnings("unchecked")
public static ContinuousDistribution getDistributionMLE
                    (String distName, double[] x, int n) {

      Class<?> distClass;
      try
      {
         distClass = Class.forName ("umontreal.ssj.probdist." + distName);
      }
      catch (ClassNotFoundException e)
      {
         try
         {
            distClass = Class.forName (distName);
         }
         catch (ClassNotFoundException ex)
         {
            throw new IllegalArgumentException ("Invalid distribution name: " 
                                                + distName);
         }
      }

      return getDistributionMLE ((Class<? extends ContinuousDistribution>)distClass, x, n);
   }
 /**
 * Uses the Java Reflection API to construct a  @ref DiscreteDistributionInt
 * object by estimating parameters of the distribution using the maximum
 * likelihood method based on the @f$n@f$ observations in table @f$x[i]@f$,
 * @f$i = 0, 1, …, n-1@f$.
 *  @param distName     the name of the distribution to instanciate
 *  @param x            the list of observations to use to evaluate
 *                      parameters
 *  @param n            the number of observations to use to evaluate
 *                      parameters
 */
@SuppressWarnings("unchecked")
public static DiscreteDistributionInt getDistributionMLE
                    (String distName, int[] x, int n) {

      Class<?> distClass;
      try
      {
         distClass = Class.forName ("umontreal.ssj.probdist." + distName);
      }
      catch (ClassNotFoundException e)
      {
         try
         {
            distClass = Class.forName (distName);
         }
         catch (ClassNotFoundException ex)
         {
            throw new IllegalArgumentException ("Invalid distribution name: " 
                                                + distName);
         }
      }

      return getDistributionMLE ((Class<? extends DiscreteDistributionInt>)distClass, x, n);
   }
 /**
 * Uses the Java Reflection API to construct a  @ref ContinuousDistribution
 * object by estimating parameters of the distribution using the maximum
 * likelihood method based on the @f$n@f$ observations in table @f$x[i]@f$,
 * @f$i = 0, 1, …, n-1@f$.
 *  @param distClass    the class of the distribution to instanciate
 *  @param x            the list of observations to use to evaluate
 *                      parameters
 *  @param n            the number of observations to use to evaluate
 *                      parameters
 */
@SuppressWarnings("unchecked")
public static <T extends ContinuousDistribution> T getDistributionMLE
                    (Class<T> distClass, double[] x, int n) {
      if (ContinuousDistribution.class.isAssignableFrom(distClass) == false)
               throw new IllegalArgumentException 
                  ("The given class is not a Probdist distribution class.");

      Method m;
      try
      {
         m = distClass.getMethod ("getInstanceFromMLE", double[].class, int.class);
      }
      catch (NoSuchMethodException e) {
         throw new IllegalArgumentException
         ("The given class does not provide the static method getInstanceFromMLE (double[],int)");
      }
      if (!Modifier.isStatic (m.getModifiers()) ||
          !distClass.isAssignableFrom (m.getReturnType()))
         throw new IllegalArgumentException
         ("The given class does not provide the static method getInstanceFromMLE (double[],int)");
      
      try
      {
         return (T)m.invoke (null, x, n);
      }
      catch (IllegalAccessException e) {
         return null;
      }
      catch (IllegalArgumentException e) {
         return null;
      }
      catch (InvocationTargetException e) {
         return null;
      }      
   }
 /**
 * Uses the Java Reflection API to construct a  @ref DiscreteDistributionInt
 * object by estimating parameters of the distribution using the maximum
 * likelihood method based on the @f$n@f$ observations in table @f$x[i]@f$,
 * @f$i = 0, 1, …, n-1@f$.
 *  @param distClass    the class of the distribution to instanciate
 *  @param x            the list of observations to use to evaluate
 *                      parameters
 *  @param n            the number of observations to use to evaluate
 *                      parameters
 */
@SuppressWarnings("unchecked")
public static <T extends DiscreteDistributionInt> T getDistributionMLE
                    (Class<T> distClass, int[] x, int n) {
      if (DiscreteDistributionInt.class.isAssignableFrom(distClass) == false)
               throw new IllegalArgumentException 
                  ("The given class is not a discrete distribution class over integers.");
      
      Method m;
      try
      {
         m = distClass.getMethod ("getInstanceFromMLE", int[].class, int.class);
      }
      catch (NoSuchMethodException e) {
         throw new IllegalArgumentException
         ("The given class does not provide the static method getInstanceFromMLE (int[],int)");
      }
      if (!Modifier.isStatic (m.getModifiers()) ||
          !distClass.isAssignableFrom (m.getReturnType()))
         throw new IllegalArgumentException
         ("The given class does not provide the static method getInstanceFromMLE (int[],int)");
      
      try
      {
         return (T)m.invoke (null, x, n);
      }
      catch (IllegalAccessException e) {
         return null;
      }
      catch (IllegalArgumentException e) {
         return null;
      }
      catch (InvocationTargetException e) {
         return null;
      }      
   }

   /**
    * Uses the Java Reflection API to construct a
    * @ref ContinuousDistribution object by executing the code contained
    * in the string `str`. This code should be a valid invocation of the
    * constructor of a  @ref ContinuousDistribution object. This method
    * throws exceptions if it cannot parse the given string and returns
    * `null` if the distribution object could not be created due to a
    * Java-specific instantiation problem.
    *  @param str          string that contains a call to the constructor
    *                      of a continuous distribution
    *  @return a continuous distribution object or `null` if it could not
    * be instantiated
    *
    *  @exception IllegalArgumentException if parsing problems occured
    * when reading `str`
    *  @exception ClassCastException if the distribution string does not
    * represent a continuous distribution
    */
   public static ContinuousDistribution getContinuousDistribution (String str) {
      return (ContinuousDistribution)getDistribution (str);
   }

   /**
    * Same as  #getContinuousDistribution, but for discrete distributions
    * over the real numbers.
    *  @param str          string that contains a call to the constructor
    *                      of a discrete distribution
    *  @return a discrete distribution object, or `null` if it could not
    * be instantiated
    *
    *  @exception IllegalArgumentException if parsing problems occured
    * when reading `str`
    *  @exception ClassCastException if the distribution string does not
    * represent a discrete distribution
    */
   public static DiscreteDistribution getDiscreteDistribution (String str) {
      return (DiscreteDistribution)getDistribution (str);
   }

   /**
    * Same as  #getContinuousDistribution, but for discrete distributions
    * over the integers.
    *  @param str          string that contains a call to the constructor
    *                      of a discrete distribution
    *  @return a discrete distribution object, or `null` if it could not
    * be instantiated
    *
    *  @exception IllegalArgumentException if parsing problems occured
    * when reading `str`
    *  @exception ClassCastException if the distribution string does not
    * represent a discrete distribution
    */
   public static DiscreteDistributionInt getDiscreteDistributionInt (String str) {
      return (DiscreteDistributionInt)getDistribution (str);
   }
}