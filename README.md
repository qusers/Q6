qsource
=======
Version control of the molecular dynamics code called Q, version 5, via github.

Q is a set of Molecular Dynamics (MD) tools tailored to specific kinds of free energy calculations, mainly:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)

## Installation
The current compilation works in Mac OSX 10.9.1 using gfortran via fink.
So, in order for it to compile you will first have to make sure that
gfortran is installed in your fink folder and also that fink has been started, i.e.:
source /sw/bin/init.sh

```bash
git clone https://github.com/qusers/qsource.git
cd qsource/src
source /sw/bin/init.sh
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

20/01/2014

In order to compile using gfortran with Alexandre Barrozo's additional subroutine for angle constraints
it's necessary to turn off the flush and iargc symbols.



