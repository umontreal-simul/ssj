/*
 * Class:        MultiFunction
 * Description:
 * Environment:  Java
 * Software:     SSJ
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       Richard Simard
 * @since        May 2013
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
package umontreal.ssj.functions;

/**
 * This interface should be implemented by classes which represent
 * *multivariate* mathematical functions. It is used to pass an arbitrary
 * function of a vector variable as argument to another function.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public interface MultiFunction {

   /**
    * Returns the value of the function evaluated at @f$X@f$.
    *  @param X            point at which the function is evaluated
    *  @return value of function at `X`
    */
   public double evaluate (double[] X);

}