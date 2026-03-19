/**
 * @package umontreal.ssj.latnetbuilder
 *
 * Provides an interface to the `LatNet Builder` software to search for good quasi-Monte Carlo point sets and sequences 
 * defined as integration lattices, polynomial lattices, and digital nets and sequences.
 * 
 * @anchor REF_latnetbuilder_overview_sec_overview
 *
 * This package provides an interface to the *LatNet Builder* software to search for good parameters for
 * quasi-Monte Carlo point sets and sequences (or HUPS) of various types 
 * (lattices, polynomial lattices, digital nets, etc.) for arbitrary dimension and cardinality, 
 * various uniformity criteria, etc. 
 * See @cite vLEC12a, @cite vLEC16a, @cite vLEC22m, and 
 * [https://github.com/umontreal-simul/latnetbuilder](https://github.com/umontreal-simul/latnetbuilder)
 * for more details on `LatNet Builder` and its ancestor `Lattice Builder`.
 * 
 * *Important:* 
 * For this package to work, `LatNet Builder` must be installed and the path to it must be set correctly
 * inside the class @ref Search.  Note that `LatNet Builder` does not run directly under *Windows*,
 * but only under Linux or in Windows in a virtual machine; see the `LatNet Builder` documentation.
 * 
 */

package umontreal.ssj.latnetbuilder;
