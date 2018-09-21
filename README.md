# SSJ

*Stochastic Simulation in Java*

SSJ is a Java library for stochastic simulation, developed under the supervision of
[Pierre L'Ecuyer](http://www-labs.iro.umontreal.ca/~lecuyer/)
in the 
[Simulation and Optimization Laboratory](http://simul.iro.umontreal.ca/),
[Department of Computer Science and Operations
Research](http://en.diro.umontreal.ca) at Université de Montréal.
It provides facilities for:

- random number and random variate generation
- stochastic process simulation
- discrete-event simulation
- computations with several types of probability distributions
- randomized quasi-Monte Carlo methods
- collecting and reporting statistics from simulations
- goodness-of-fit tests
- and much more.

Starting from version 3.1.0, SSJ is released under the
[Apache 2.0 License](http://www.apache.org/licenses/LICENSE-2.0),
and the package names have changed from `umontreal.iro.lecuyer.*` to `umontreal.ssj.*`.


## Documentation and tutorial

The [SSJ User's Guide](http://umontreal-simul.github.io/ssj/docs/master)
includes:
- the [API documentation](http://umontreal-simul.github.io/ssj/docs/master/namespaces.html); and
- [a tutorial with documented examples](http://umontreal-simul.github.io/ssj/docs/master/examples.html).


## Installation

You can install SSJ either by [adding it as a dependency](#using-maven)
for your Maven- or Gradle-based project, by downloading a [binary
release](#binary-releases) or by [compiling it from
scratch](#compiling-the-source-code).

SSJ is compatible with Java SE8 (JDK 8) and later versions of Java. The latest
Java JDK is available at
[Oracle](http://www.oracle.com/technetwork/java/javase/downloads/index.html)
with installation instructions.  It should be installed *before* installing SSJ.

### Using Maven

**(Simplest approach)**

SSJ packages are hosted
[on Bintray](https://bintray.com/umontreal-simul/maven/ssj/_latestVersion)
and on
[Maven Central](https://repo1.maven.org/maven2/ca/umontreal/iro/simul/ssj/).
If your Java project uses [Maven](http://maven.apache.org/) or
[Gradle](http://gradle.org/), **all you need to do is add
`ca.umontreal.iro.simul:ssj:+` to the Maven dependencies of your project**,
then you can start working on your SSJ-based project right-away.


##### IDE integration

Several integrated development environments (IDEs) such as 
[Eclipse](http://www.eclipse.org/),
[NetBeans](http://netbeans.org/),
[IntelliJ IDEA](http://www.jetbrains.com/idea/), for example, support Maven.  
If you work in one of these IDEs, simply create your project as a *Maven project* instead
of a Java project, then add SSJ to the Maven dependencies of the
project (refer to your IDE documentation), with the following parameters:

- **Group Id**: `ca.umontreal.iro.simul`;
- **Artifact Id**: `ssj`;
- **Version**: any [valid release number for SSJ](http://github.com/umontreal-simul/ssj/releases),
  e.g., `3.3.0`,

and you are ready to go!


### Binary releases

For those who want to download the binaries and install them manually, 
we provide below some general instructions for configuring a project to use SSJ.
Less experienced users can find
[more detailed instructions on the SSJ page](http://simul.iro.umontreal.ca/ssj/indexe.html).


#### Download a binary archive

Pre-compiled binaries are available as archives on the
[releases page](http://github.com/umontreal-simul/ssj/releases).
They include:

- the main SSJ [JAR file](http://en.wikipedia.org/wiki/JAR_%28file_format%29)
  (under `ssj/lib`);
- JAR files for the Java libraries used by SSJ (the dependencies)
  (under `ssj/lib`);
- dependencies and the [JNI](http://en.wikipedia.org/wiki/Java_Native_Interface)
  shared libraries (under the `ssj/lib` directory);
- the user's guide (under `ssj/doc/html`); and
- example source files (under `ssj/doc/examples`).

You can download the latest archive and extract the files in a location of your choice.


#### Set the Java class path

You add to the Java class path every JAR file found under the `ssj/lib`
directory of the binary archive.

##### On the command line

If you use Java from the command line, add the full path of every JAR file
under `ssj/lib` to the `CLASSPATH` environment variable, separated with `:`
under Linux or MacOS, or by `;` under Windows.
Means of doing this depends on the system you are using.
For example, under Linux with a Bash-compatible shell, one could use something
like:

```sh
for f in /full/path/to/ssj/lib/*.jar; do
    CLASSPATH=$f:$CLASSPATH
done
export CLASSPATH
```


##### In Eclipse

In Eclipse, under *Window ‣ Preferences ‣ Java ‣ Build Path ‣ User Libraries*,
click *New…*.  Set the name to `SSJ` and click *OK*.  Click *Add External
JARs…* navigate to the `ssj/lib` folder of the extracted binary archive, select
all JAR files, and click *OK*.  You can now add the SSJ library you have
created from any project, by right-clicking on your project name in the
*Package Explorer*, by selecting *Build Path ‣ Add Libraries… ‣ User Library*
under your project tree and by choosing SSJ.


##### In NetBeans

In NetBeans, under *Tools ‣ Libraries*, press *New Library…*.  Set the name to
`SSJ` and click *OK*.  Click *Add JAR/Folder…*, navigate to the `ssj/lib`
folder of the extracted binary archive, select all JAR files, and click *Add
JAR/Folder*.  You can now add the SSJ library you have created from any
project, by right-clicking on *Libraries* under your project tree in the
*Projects* tab and by choosing SSJ.



## Compiling the source code

**You do not need to compile the source code to use SSJ if you have already installed it
[using Maven](#using-maven) or a [binary release](#binary-releases).**
But in case you want to change the source for some reason, here is how you can recompile.

The SSJ library uses [Gradle](http://gradle.org/) as its build system.
You do not need to download it, since the Gradle wrapper executable program is
provided with the source code as the `gradlew` file for Linux and MacOS
platforms, and as `gradlew.bat` for Windows platforms.
The `build.gradle` and `gradle.properties` files at the root of the source tree
contain the configuration for Gradle.
In the instructions below, **Windows users** should replace instances of
`./gradlew` with `gradlew.bat`.

SSJ and the current Gradle version work with Java SE (or JDK) version 8 or later.


### Using Gradle

The general syntax for Gradle is `./gradlew <task>` where `<task>` is the name
of a Gradle task (run `./gradlew tasks` to obtain a list of available tasks).

On the command line, from the root of the source tree, type:

- `./gradlew check` to build and test the library;
- `./gradlew examples` (optionally) to run additional examples;
- `./gradlew distZip` or `./gradlew distTar` to create an binary archive of the
  SSJ library, including the SSJ JAR file and its dependencies;
- `./gradlew --gui` to launch the Gradle graphical user interface and choose
  from more options.

All files generated during the build process are placed under the `build`
subdirectory.
The generated binary archives can be found under `build/distributions`.


### Building the documentation

**(Optional)**

The SSJ library uses [Doxygen](http://www.stack.nl/~dimitri/doxygen/) as its
documentation system.
If Doxygen is available on your system, you can tell Gradle to build the
documentation by adding the following line in `gradle.properties`:

    buildDocs

Then, run Gradle [as explained above](#using-gradle).
You may want (or need) to change some Doxygen environments variables in the file Doxyfile.


### JNI classes

**(Optional)**

The classes
[UnuranContinuous](http://umontreal-simul.github.io/ssj/docs/master/html/classumontreal_1_1ssj_1_1randvar_1_1UnuranContinuous.html),
[UnuranDiscreteInt](http://umontreal-simul.github.io/ssj/docs/master/html/classumontreal_1_1ssj_1_1randvar_1_1UnuranDiscreteInt.html),
[UnuranEmpirical](http://umontreal-simul.github.io/ssj/docs/master/html/classumontreal_1_1ssj_1_1randvar_1_1UnuranEmpirical.html) and
[GlobalCPUTimeChrono](http://umontreal-simul.github.io/ssj/docs/master/html/classumontreal_1_1ssj_1_1util_1_1GlobalCPUTimeChrono.html)
make use of native libraries through the
[Java Native Interface (JNI)](http://en.wikipedia.org/wiki/Java_Native_Interface).
These libraries must be compiled with a C compiler (known to work with GCC).

Note that if you want to build and use the UNU.RAN interface provided with SSJ,
you must first install [UNU.RAN](http://statistik.wu-wien.ac.at/unuran/).

To tell Gradle to build the JNI libraries, add the following lines in
`gradle.properties`:

    ssjutil.jni.build
    randvar.jni.build
    unuran.prefix = "/path/to/unuran"

Make sure to replace `/path/to/unuran` in the above with the installation
prefix of UNU.RAN, i.e., the path under which `include/unuran.h` can be found.
Then, run Gradle [as explained above](#using-gradle).


### Cross-compiling

**(For experts only!)**

Under Linux with [GCC](http://gcc.gnu.org/) and [MinGW](http://www.mingw.org/),
you can cross-compile the JNI libraries for both Linux and Windows.
This is how we generate the binary archives.

To enable cross-compilation with Gradle, add the following lines in
`gradle.properties`:

    crossCompile
    unuran.prefix.linux64 = "/path/to/unuran"
    unuran.prefix.win32   = "/path/to/unuran-mingw32"

Make sure to replace `/path/to/unuran` and `/path/to/unuran-mingw32` with the
installation prefixes of UNU.RAN compiled for 64-bit Linux and for 32-bit
Windows, respectively.


## Dependencies

SSJ depends on the following libraries.
**You do not need to download them manually if you're using SSJ with Maven or
from a binary release.**

If you intend to [compile the source code](#compiling-the-source-code) of SSJ,
Gradle will take care of downloading the Java dependencies for you.
Optionally, if you want to use the UNU.RAN interface in SSJ, you need to
install the [UNU.RAN](http://statistik.wu-wien.ac.at/unuran/) before compiling
the associated JNI shared library in SSJ.

##### [Colt](http://dst.lbl.gov/ACSSoftware/colt/)  
The Colt library is used by a few SSJ classes.  The library, its source code
and documentation, can be downloaded for free from its home page.  The
`colt.jar` archive is already included in the SSJ distribution and it must be
in the CLASSPATH environment variable.

##### [Nonlinear Optimization Java Package](http://www1.fpl.fs.fed.us/optimization.html) by [Steve Verrill](http://www1.fpl.fs.fed.us/steve.html)  
The optimization package of Steve Verrill includes Java translations of the
MINPACK nonlinear least squares routines as well as UNCMIN routines for
unconstrained optimization.  They were translated from FORTRAN to Java by
Steve Verrill and are in the public domain.  They are included in the SSJ
distribution as the optimization.jar archive.  It is used only in the
`probdist` package to compute maximum likelihood estimators.

##### [JFreeChart](http://www.jfree.org/jfreechart/)
The JFreeChart library is used by the SSJ package charts to draw curves, histograms
and different kinds of plots.  JFreeChart is copyrighted under the
[GNU LGPL License](http://www.gnu.org/licenses/lgpl.html).
It is included in the SSJ distribution as the `jfreechart-X.Y.Z.jar`
and the `jcommon-X.Y.Z.jar`, where `X.Y.Z` represents a version number.

##### [UNU.RAN](http://statistik.wu-wien.ac.at/unuran/) (optional)  
The UNURAN library is used by the classes `UnuranContinuous`,
`UnuranDiscrete` and `UnuranEmpirical` in the package called `randvar`.
Downloading, compiling and installing UNURAN is optional.  It is required
only if SSJ must be rebuilt.  However, the UNURAN documentation is required
to use the SSJ UNURAN interface efficiently.
