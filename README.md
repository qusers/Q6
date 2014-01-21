qsource
=======
Version control of the molecular dynamics code called Q, version 5, via github.

## Installation



```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
make all
make clean
```

Q is a set of Molecular Dynamics (MD) tools tailored to specific kinds of free energy calculations, mainly:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)


The current compilation works in Mac OSX 10.9.1 using gfortran via fink.
So, in order for it to compile you will first have to make sure that
gfortran is installed in your fink folder and also that fink has been started, i.e.:
source /sw/bin/init.sh

After that you have to cd to the src folder andn just type:

make all

This will take care of redirecting the binaries and object files to folders for code tidyness.
It is also good to do:

make clean

After compilation as that will clean the src folder.


NOTES:
=========

20/01/2014

In order to compile, using gfortran,  with Alex's additional subroutine for angle constraints it's necessary to
turn of the symbols flush and iargc, which he has turned on as he's most likely compiling with ifortran.
