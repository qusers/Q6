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
make mpi COMP=gcc
```

To compile in Mac OSX 10.9.2 you can use native gfortran binaries which you can download from  (https://gcc.gnu.org/wiki/GFortranBinaries#MacOS) or you can also compile using the GCC (Gnu Compiler Collection) distributions available to fink, macports or homebrew. In order to compile in a MAC you should call the fink environment first, usually with:

```bash
source /sw/bin/init.sh
```

alternatively, you can install [homebrew](http://brew.sh/) and use (confirmed for Mac OSX 10.10):
```bash
brew install gcc
```


Then you proceed as in linux
```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
make all COMP=osx
```

This will take care of redirecting the binaries and object files to standard bin and obj folders for code tidyness.

After this you have to add the program to your system path by modifying your shell initiation script, that is, if your shell is bash, you can add the following lines to your .bashrc file using a text editor:

```bash
SOFT=/Users/johndoe/software
export QDIR=$SOFT/qsource
export PATH=$QDIR/bin:$PATH  
```
Where $SOFT will be the place where your software folder is located at, e.g. /Users/johndoe/software

Once the q binaries are declared in your path you should be able to call all q binaries from your terminal.
To test that the path to your compiled Q binaries has been correctly assigned you can issue the following commands in the terminal:
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


