/*
 * Class:        ClassFinder
 * Description:  Convert a simple class name to a fully qualified class object
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

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Utility class used to convert a simple class name to a fully qualified
 * class object. The  Class class can be used to obtain information about a
 * class (its name, its fields, methods, constructors, etc.), and to
 * construct objects, even if the exact class is known at runtime only. It
 * provides a  java.lang.Class.forName static method converting a string to a
 * Class, but the given string must be a fully qualified name.
 *
 * Sometimes, configuration files may need to contain Java class names. After
 * they are extracted from the file, these class names are given to
 * java.lang.Class.forName to be converted into  Class objects.
 * Unfortunately, only fully qualified class names will be accepted as input,
 * which clutters configuration files, especially if long package names are
 * used. This class permits the definition of a set of import declarations in
 * a way similar to the Java Language Specification @cite iGOS00a&thinsp;. It
 * provides methods to convert a simple class name to a  Class object and to
 * generate a simple name from a  Class object, based on the import rules.
 *
 * The first step for using a class finder is to construct an instance of
 * this class. Then, one needs to retrieve the initially empty list of import
 * declarations by using  #getImports, and update it with the actual import
 * declarations. Then, the method  #findClass can find a class using the
 * import declarations. For example, the following code retrieves the class
 * object for the `List` class in package `java.util`
 *
 * @code
 *
 *    ClassFinder cf = new ClassFinder();
 *    cf.getImports().add ("java.util.*");
 *    Class<?> listClass = cf.findClass ("List");
 *
 * @endcode
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public class ClassFinder implements Cloneable, java.io.Serializable {
   private static final long serialVersionUID = -4847630831331065792L;

/**
 * Contains the saved import lists. Each element of this list is a nested
 * List containing  String ’s, each string containing the fully qualified
 * name of an imported package or class.
 *  @serial
 */
private List<List<String>> imports = new LinkedList<List<String>>();

   /**
    * Constructs a new class finder with an empty list of import
    * declarations.
    */
   public ClassFinder() {
      List<String> imp = new ArrayList<String>();
      imports.add (imp);
   }

   /**
    * Returns the current list of import declarations. This list may
    * contain only  String ’s of the form `java.class.name` or
    * `java.package.name.*`.
    *  @return the current list of import declarations.
    */
   public List<String> getImports() {
      return imports.get (imports.size() - 1);
   }

   /**
    * Saves the current import list on the import stack. This method makes
    * a copy of the list returned by  #getImports and puts it on top of a
    * stack to be restored later by  #restoreImports.
    */
   public void saveImports() {
      List<String> imp = getImports();
      List<String> impBack = new ArrayList<String> (imp);
      imports.add (impBack);
   }

   /**
    * Restores the list of import declarations. This method removes the
    * last list of import declarations from the stack. If the stack
    * contains only one list, this list is cleared.
    */
   public void restoreImports() {
      if (imports.size() == 1)
         getImports().clear();
      else
         imports.remove (imports.size() - 1);
   }

   /**
    * Tries to find the class corresponding to the simple name `name`. The
    * method first considers the argument as a fully qualified class name
    * and calls  java.lang.Class.forName `(name)`. If the class cannot be
    * found, it considers the argument as a simple name. A simple name
    * refers to a class without specifying the package declaring it. To
    * convert simple names to qualified names, the method iterates through
    * all the strings in the list returned by  #getImports, applying the
    * same rules as a Java compiler to resolve the class name. However, if
    * an imported package or class does not exist, it will be ignored
    * whereas the compiler would stop with an error.
    *
    * For the class with simple name `name` to be loaded, it must be
    * imported explicitly (single-type import) or one of the imported
    * packages must contain it (type import on-demand). If the class with
    * name `name` is imported explicitly, this import declaration has
    * precedence over any imported packages. If several import declaration
    * match the given simple name, e.g., if several fully qualified names
    * with the same simple name are imported, or if a class with simple
    * name `name` exists in several packages, a
    * @ref NameConflictException is thrown.
    *  @param name         the simple name of the class.
    *  @return a reference to the class being loaded.
    *
    *  @exception ClassNotFoundException if the class cannot be loaded.
    *  @exception NameConflictException if a name conflict occurred.
    */
   public Class<?> findClass (String name) throws
      ClassNotFoundException, NameConflictException {
      // Try to consider the name as a fully qualified class name
      try {
         return Class.forName (name);
      }
      catch (ClassNotFoundException cnfe) {}

      List<String> imports = getImports();
      Class<?> candidate = null;
      String candidateImportString = "";
      boolean candidateImportOnDemand = false;

      // Determines the name of the outermost class
      // if name corresponds to the simple name of a nested class, e.g., for A.B.C
      // outerName will contain A.
      int idxOut = name.indexOf ('.');
      String outerName;
      if (idxOut == -1)
         outerName = name;
      else
         outerName = name.substring (0, idxOut);
      for (String importString : imports) {
         // For each import declaration, we try to load
         // a class and store the result in cl.
         // When cl is not null, the Class object is
         // compared with the best candidate candidate.
         Class<?> cl = null;
         boolean onDemand = false;
         if (!importString.endsWith (".*")) {
            // Single-type import declaration
            // We must ensure that this will correspond to
            // the desired class name. For example, if we
            // search List and the import java.util.ArrayList
            // is found in importString, the ArrayList
            // class must not be returned.
            if (importString.endsWith ("." + outerName)) {
               // The name of outer class was found in importString and
               // a period is found left to it.
               // So try to load this class.
               // Simple class names have precedence over
               // on-demand names.

               // Replace, in importString, the name of
               // the outer class with the true class name
               // we want to load.  If the simple name
               // does not refer to a nested class, this
               // has no effect.
               String cn = importString.substring
                  (0, importString.length() - outerName.length()) + name;
               try {
                  cl = Class.forName (cn);
               }
               catch (ClassNotFoundException cnfe) {}
            }
         }
         else {
            // Type import on demand declaration
            try {
               // Replace the * with name and
               // try to load the class.
               // If that succeeds, our candidate cl
               // is onDemand.
               cl = Class.forName
                  (importString.substring (0, importString.length() - 1) + name);
               onDemand = true;
            }
            catch (ClassNotFoundException cnfe) {}
         }
         if (cl != null) {
            // Something was loaded
            if (candidate == null ||
                (candidateImportOnDemand && !onDemand)) {
               // We had no candidate or the candidate was imported
               // on-demand while this one is a single-type import.
               candidate = cl;
               candidateImportString = importString;
               candidateImportOnDemand = onDemand;
            }
            else if (candidate != cl)
               throw new NameConflictException
                  (this, name,
                   "simple class name " + name +
                   " matches " + candidate.getName() +
                   " (import string " + candidateImportString + ") or " +
                   cl.getName() + " (import string " + importString + ")");
         }
      }
      if (candidate == null)
         throw new ClassNotFoundException
            ("Cannot find the class with name " + name);
      return candidate;
   }

   /**
    * Returns the simple name of the class `cls` that can be used when the
    * imports contained in this class finder are used. For example, if
    * `java.lang.String.class` is given to this method, `String` is
    * returned if `java.lang.*` is among the import declarations.
    *
    * Note: this method does not try to find name conflicts. This
    * operation is performed by  #findClass(String) only. For example, if
    * the list of imported declarations contains `foo.bar.*` and
    * `test.Foo`, and the simple name for `test.Foo` is queried, the
    * method returns `Foo` even if the package `foo.bar` contains a `Foo`
    * class.
    *  @param cls          the class for which the simple name is queried.
    *  @return the simple class name.
    */
   public String getSimpleName (Class<?> cls) {
      if (cls.isArray())
         return getSimpleName (cls.getComponentType()) + "[]";
      if (cls.isPrimitive())
         return cls.getName();
      Class<?> outer = cls;
      while (outer.getDeclaringClass() != null)
         outer = outer.getDeclaringClass();
      boolean needsFullyQualified = true;
      for (String importString : getImports()) {
         if (importString.equals (outer.getName()))
             // A single-type import is given, can return an unqualified name.
            needsFullyQualified = false;
         else if (importString.endsWith (".*")) {
            // Remove the .* at the end of the import string to get a package name.
            String pack = importString.substring (0, importString.length() - 2);
            // Compare this package name.
            if (pack.equals (cls.getPackage().getName()))
               needsFullyQualified = false;
         }
      }
      if (needsFullyQualified)
         return cls.getName();
      else {
         String name = cls.getName();
         String pack = cls.getPackage().getName();
         if (!name.startsWith (pack))
            throw new IllegalStateException
               ("The class name " + name +
                " does not contain the package name " + pack);

         // Removes the package and the . from the fully qualified class name.
         return name.substring (pack.length() + 1);
      }         
   }

   /**
    * Clones this class finder, and copies its lists of import
    * declarations.
    */
   public ClassFinder clone() {
      ClassFinder cf;
      try {
         cf = (ClassFinder)super.clone();
      }
      catch (CloneNotSupportedException cne) {
         throw new InternalError
            ("CloneNotSupported thrown for a class implementing Cloneable");
      }
      cf.imports = new LinkedList<List<String>>();
      for (List<String> imp : imports) {
         List<String> impCpy = new ArrayList<String> (imp);
         cf.imports.add (impCpy);
      }
      return cf;
   }
}