qsource
=======
Version control of the molecular dynamics code called Q, version 5, via github.

Q is a set of Molecular Dynamics (MD) tools tailored to specific kinds of free energy calculations, mainly:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)

## Installation
The current makefiles make it relatively easy to compile the code in a linux or mac OSX environment.
You have to make sure first that gfortran is installed and if you want to compile the MPI version
you also have to make sure that openmpi is installed.

To install in a linux environment you have to move to the src/ folder and rename the correponding
makefile.
```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
cp makefile.linux makefile
make all
make clean
```

Or for MPI compilation:
```bash
make mpi
```

As the included makefiles use gfortran and/or mpif90 one needs to use fink, macports or homebrew in Mac's
OSX 10.9.1. We've tested with fink, but it should compile without major problems using the others.
In order to compile in a MAC you should call the fink environment first, usually with:

```bash
source /sw/bin/init.sh
```

Then you proceed as in linux
```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
cp makefile.osx makefile
make all
make clean
```

This will take care of redirecting the binaries and object files to folders for code tidyness.

After this you have to add the program to your system path by modifying your shell initiation script, that is,
if your shell is bash you can add the following line to your .bashrc file:

```bash
export Q5MD=$SOFT/qsource
export PATH=$Q5MD/bin:$PATH  
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

29/01/2014
Changed fortran source to Paul's current version.

20/01/2014

In order to compile using gfortran with Alexandre Barrozo's additional subroutine for angle constraints
it was necessary to turn off the flush and iargc symbols.

