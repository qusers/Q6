qsource
=======
Version control of the molecular dynamics code called Q, version 5.06, via github.

Q is a set of Molecular Dynamics (MD) tools tailored to specific kinds of free energy calculations, mainly:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)

## Installation
The current makefiles make it relatively easy to compile the code in a linux or mac OSX environment.

You have to make sure first that gfortran is installed and if you want to compile the MPI version
you also have to make sure that openmpi is installed.

If you're using gfortran make sure that you have version 4.8 or later. This applies both to compilation in mac and linux.
To make sure that you have version 4.8 or later use:

```bash
gfortran --version
```

Right before issuing the "make" command.

To install in a linux environment you have to move to the src/ folder where the source code and the makefile are located at. To get information on how to use the makefile just type make in your terminal:
```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
make
```
Following the instructions for compilation using the gfortran compiler and the non-mpi version of the code would then take the form:

```bash
make all COMP=gcc
```

Or for MPI compilation (after loading a proper MPI compiler library):
```bash
cp qdyn.F90_ifort_signals qdyn.f90
make mpi COMP=gcc
```

As the included makefiles use gfortran and/or mpif90 one needs to use fink, macports or homebrew in Mac's OSX 10.9.2. We've tested with fink, but it should compile without major problems using the others.
In order to compile in a MAC you should call the fink environment first, usually with:

```bash
source /sw/bin/init.sh
```

Then you proceed as in linux
```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
make all COMP=osx
```

This will take care of redirecting the binaries and object files to standard bin and obj folders for code tidyness.

After this you have to add the program to your system path by modifying your shell initiation script, that is, if your shell is bash, you can add the following line to your .bashrc file:

```bash
export QDIR=$SOFT/qsource
export PATH=$QDIR/bin:$PATH  
```
Where $SOFT will be the place where your software folder is located at, e.g. /Users/johndoe/software

Once the q binaries are declared in your path you should be able to call all q binaries from your terminal.

```bash
source .bashrc
echo $path | grep qsource

/Users/johndoe/software/qsource

qprep5

###############################################################################
Welcome to Qprep version 5.04

Qprep> 
```



NOTES:
=========

25/04/2014

A new makefile is ready allowing easy upgrades/changes on compiler options by just adding additional else if statements. Work on updating headers and version numbering has been added to the code. Author information added to the main makefile. The script for sending tests has been fixed to correctly account for the number of cores per cluster node.


10/02/2014

If one compiles with older versions of gfortran, say, 4.6, a compilation error will come up when compiling the qcalc program at the rdf file.
This issue is not present in the ifortran compilation at all as the fortran code seemingly makes good use of intel fortran options.


07/02/2014

Checked that the code compiles in Mac OSX Mavericks natively using Xcode developer tools and compiled recent versions of gcc and gfortran from:

http://goo.gl/pk5nq5


29/01/2014

Changed fortran source to Paul's current version.


20/01/2014

In order to compile using gfortran with Alexandre Barrozo's additional subroutine for angle constraints
it was necessary to turn off the flush and iargc symbols.


