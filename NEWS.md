# SSJ release highlights
---


## SSJ 3.3.0

This release contains new packages, several new classes, an improved tutorial and improved documentation in general.

- new packages: `mcqmctools`, `latnetbuilder`, `stat.density`, `discrepancy`.
- several methods have been added in many places.
- the main API page (introduction), the tutorial (with examples), and the documentation of several packages have been updated.

## SSJ 3.2.1

This release contains bug fixes.

- package `functionfit`: Fix bug when creating an approximate BSpline.
- package `util`: bisection and brentDekker methods will now check the bounds of the interval first.
- package `eventlist`: Removed free node stack and removed synchronization on the free node stack. This increases the performance when executing many simulations in multi-thread program.

## SSJ 3.2.0

- Add new packages `markovchainrqmc` and `util.sort`.
- Add new packages `stat.list.lincv` and `stat.matrix`.
- Add new classes to package `stat` and `stat.list`.
- Add new classes to package `hups`.
- Add Rijndael's algorithm to package `rng`.
- Add new classes to package `stochprocess`.

## SSJ 3.0.0-rc1

This version of SSJ is not for production.

We have just migrated the source code from a custom documentation system which generated Java files from LaTeX input files. Now the Java source code is readily editable. The documentation is in now Doxygen format, which supports mathematical formulas and BibTeX citations.

Users of SSJ 2.6 who want to upgrade their code for SSJ 3 need to replace umontreal.iro.lecuyer with umontreal.ssj in the package imports. For example,

    import umontreal.iro.lecuyer.rng.MRG32k3a; // using SSJ 2

becomes

    import umontreal.ssj.rng.MRG32k3a; // using SSJ 3


